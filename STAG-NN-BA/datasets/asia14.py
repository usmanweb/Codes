import glob
import os
import os.path as osp

import torch
from torch import LongTensor
from torch.utils.data import random_split
from torch_geometric.data import Data, InMemoryDataset
from torch_geometric.data import download_url
from torchvision import datasets
from tqdm import tqdm


def load_asia14(data_dir, split=0.2, transform=None, seed=None):
    if seed:
        torch.manual_seed(seed)

    root = str(os.path.join(data_dir, 'Training Set'))
    asia14 = datasets.ImageFolder(root=root, transform=transform)
    lengths = [round(len(asia14) * (1-split)), round(len(asia14) * split)]
    train, test = random_split(asia14, lengths)
    return train, test


class ASIA14(InMemoryDataset):
    url = 'https://drive.google.com/u/0/uc?id=1-IdZu70SqCMLa5kBGNnuZJzjIQAKoYqQ&export=download'

    def __init__(self, root, train=True, split=0.2, transform=None, pre_transform=None,
                 pre_filter=None, seed=None):
        self.split = split
        self.seed = seed
        super().__init__(root, transform, pre_transform, pre_filter)
        path = self.processed_paths[0] if train else self.processed_paths[1]
        self.data, self.slices = torch.load(path)

    @property
    def raw_dir(self) -> str:
        return osp.join(self.root)

    @property
    def raw_paths(self):
        return glob.glob(f'{self.raw_dir}/*/*/*[(.png)(.jpg)]')

    @property
    def processed_file_names(self):
        return ['train_data.pt', 'test_data.pt']

    def download(self):
        # Download to `self.raw_dir`.
        download_url(self.url, self.raw_dir)

    def process(self):
        train, test = load_asia14(self.raw_dir, split=self.split, transform=self.pre_transform, seed=self.seed)

        for split, subset in zip(['train', 'test'], [train, test]):
            data_list = [Data(**data.stores[0], y=LongTensor([y])) for data, y in tqdm(subset, total=len(subset))]
            data, slices = self.collate(data_list)
            torch.save((data, slices), osp.join(self.processed_dir, f'{split}_data.pt'))
