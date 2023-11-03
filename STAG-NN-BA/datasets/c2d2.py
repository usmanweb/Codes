import glob
import os.path as osp
import pickle

import torch
from torch import LongTensor
from torch_geometric.data import Data, InMemoryDataset, Dataset
from torch_geometric.data import download_url
from tqdm import tqdm


class C2D2Dataset(Dataset):

    def __init__(self, root, train=True, transform=None):
        super(C2D2Dataset, self).__init__()
        filename = 'train' if train else 'test'

        with open(f'{root}/{filename}.pkl', 'rb') as handle:
            dx, dy = pickle.load(handle)

        self.data = dx
        self.targets = torch.LongTensor(dy)
        self.transform = transform

    def __getitem__(self, index):
        x = self.data[index]
        y = self.targets[index]

        if self.transform:
            x = self.transform(x)

        return x, y

    def __len__(self):
        return len(self.data)


def load_c2d2(data_dir, transform=None, seed=None):
    if seed:
        torch.manual_seed(seed)

    root = str(data_dir)

    test = C2D2Dataset(root=root, train=False, transform=transform)
    train = C2D2Dataset(root=root, train=True, transform=transform)
    # lengths = [round(len(asia14) * (1-split)), round(len(asia14) * split)]
    # train, test = random_split(asia14, lengths)
    return train, test


class C2D2(InMemoryDataset):
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
        train, test = load_c2d2(self.raw_dir, transform=self.pre_transform, seed=self.seed)

        for split, subset in zip(['test', 'train'], [test, train]):
            data_list = [Data(**data.stores[0], y=LongTensor([y])) for data, y in tqdm(subset, total=len(subset))]
            data, slices = self.collate(data_list)
            torch.save((data, slices), osp.join(self.processed_dir, f'{split}_data.pt'))
