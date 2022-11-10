import pickle


def load_pickle(path):
    return pickle.load(open(path, "rb"))


class Transformers:
    def __init__(self):
        
#         self._scaler = load_pickle("./content/scaler.pkl")
        self._pca = load_pickle("./content/pca.pkl")
        self._umap = load_pickle("./content/umap.pkl")

    @property
    def scaler(self):
        return self._scaler

    @property
    def pca(self):
        return self._pca

    @property
    def umap(self):
        return self._umap