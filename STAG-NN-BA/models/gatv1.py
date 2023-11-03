import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Module
from torch_geometric.nn import GATConv, global_add_pool


class GATV1(Module):
    def __init__(self, num_features, num_classes, layer_sizes, attn_heads):
        super(GATV1, self).__init__()  # Init parent
        self.GAT_layer_sizes = [num_features] + layer_sizes
        self.layer_heads = [1] + attn_heads
        self.MLP_layer_sizes = [self.GAT_layer_sizes[-1] * attn_heads[-1], 32, num_classes]
        self.MLP_acts = [F.relu, lambda x: x]

        self.GAT_layers = nn.ModuleList(
            [
                GATConv(dim_in * heads_in, dim_out, heads=heads_out)
                for dim_in, dim_out, heads_in, heads_out in zip(
                    self.GAT_layer_sizes[:-1],
                    self.GAT_layer_sizes[1:],
                    self.layer_heads[:-1],
                    self.layer_heads[1:],
                )
            ]
        )
        self.MLP_layers = nn.ModuleList(
            [nn.Linear(dim_in, dim_out) for dim_in, dim_out in zip(self.MLP_layer_sizes[:-1], self.MLP_layer_sizes[1:])]
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        for gat_layer in self.GAT_layers:
            x = gat_layer(x=x, edge_index=edge_index)

        x = global_add_pool(x=x, batch=batch)

        for layer, activation in zip(self.MLP_layers, self.MLP_acts):
            x = activation(layer(x))

        return x


if __name__ == "__main__":
    model = GATV1(num_features=3, num_classes=10, layer_sizes=[32, 64, 64], attn_heads=[3, 3, 3])