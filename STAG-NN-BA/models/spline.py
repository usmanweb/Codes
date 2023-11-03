import torch
import torch.nn.functional as F
from torch_geometric.nn import SplineConv, global_mean_pool, graclus, max_pool, max_pool_x
from torch_geometric.utils import normalized_cut


def normalized_cut_2d(edge_index, pos):
    row, col = edge_index
    edge_attr = torch.norm(pos[row] - pos[col], p=2, dim=1)
    return normalized_cut(edge_index, edge_attr, num_nodes=pos.size(0))


class SplineNet(torch.nn.Module):
    def __init__(self, num_features, num_classes, transform=None):
        super().__init__()
        self.transform = transform
        self.conv1 = SplineConv(num_features, 32, dim=2, kernel_size=5)
        self.conv2 = SplineConv(32, 64, dim=2, kernel_size=5)
        self.fc1 = torch.nn.Linear(64, 128)
        self.fc2 = torch.nn.Linear(128, num_classes)

    def forward(self, data):
        data.x = F.elu(self.conv1(data.x, data.edge_index, data.edge_attr))
        weight = normalized_cut_2d(data.edge_index, data.pos)
        cluster = graclus(data.edge_index, weight, data.x.size(0))
        data.edge_attr = None
        data = max_pool(cluster, data, transform=self.transform)

        data.x = F.elu(self.conv2(data.x, data.edge_index, data.edge_attr))
        weight = normalized_cut_2d(data.edge_index, data.pos)
        cluster = graclus(data.edge_index, weight, data.x.size(0))
        x, batch = max_pool_x(cluster, data.x, data.batch)

        x = global_mean_pool(x, batch)
        x = F.elu(self.fc1(x))
        x = F.dropout(x, training=self.training)
        return F.log_softmax(self.fc2(x), dim=1)


def train_spline(model, data_loader, optimizer, device, epoch):
    model.train()

    if epoch == 16:
        for param_group in optimizer.param_groups:
            param_group['lr'] = 0.001

    if epoch == 26:
        for param_group in optimizer.param_groups:
            param_group['lr'] = 0.0001

    loss = 0

    for data in data_loader:
        data = data.to(device)
        optimizer.zero_grad()
        predictions = model(data)
        loss_t = F.nll_loss(predictions, data.y)
        loss += loss_t.cpu().detach().numpy()
        loss_t.backward()
        optimizer.step()

    return loss
