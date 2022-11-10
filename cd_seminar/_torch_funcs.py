
# -- import packages: --------------------------------------------------------------------
import torch
import numpy as np

# -- define functions: -------------------------------------------------------------------
def autodevice(idx: int = 0):

    """
    Return torch device.

    Parameters:
    -----------
    idx
        type: int
        default: 0

    Returns:
    --------
    torch.device

    Notes:
    ------
    (1) Looks for CUDA then mps backends and finally defaults to cpu.
    """

    if torch.cuda.is_available():
        return torch.device("cuda:{}".format(torch.cuda.current_device()))
    if torch.backends.mps.is_available():
        return torch.device("mps:{}".format(idx))
    return torch.device("cpu")




# -- define model: ---------
from collections import OrderedDict
import torch
import matplotlib.pyplot as plt
import vinplots

class GeneratorModel(torch.nn.Module):
    def __init__(self, ckpt_path=None, plot=False):
        super().__init__()
        self.net = torch.nn.Sequential(
            OrderedDict(
                [
                    ("linear1", torch.nn.Linear(50, 400)),
                    ("activation1", torch.nn.LeakyReLU()),
                    ("linear2", torch.nn.Linear(400, 400)),
                    ("activation2", torch.nn.LeakyReLU()),
                    ("linear", torch.nn.Linear(400, 1, bias=False)),
                ]
            )
        )
        
        y = self.state_dict()["net.linear1.weight"].numpy().flatten()
        w_col="lightgrey"
        if ckpt_path:
            ckpt = torch.load(ckpt_path, map_location=torch.device("cpu"))[
                "model_state_dict"
            ]
            self.load_state_dict(ckpt)
            y = self.state_dict()["net.linear1.weight"].numpy().flatten()
            w_col = "crimson"

        if plot:
            fig, axes = vinplots.quick_plot(nplots=1, figsize=0.5)
            b0 = axes[0].hist(y, bins=50, zorder=1, color=w_col, alpha=1)
            axes[0].set_title("Layer 1 Weights")
            
        print(self)
        
    def _step(self, x, dt, z):
        sqrtdt = np.sqrt(dt)
        return x + self._drift(x) * dt + z * sqrtdt

    def _pot(self, x):
        return self.net(x)

    def _drift(self, x):
        x_ = x.requires_grad_()
        pot = self._pot(x_)

        drift = torch.autograd.grad(pot, x_, torch.ones_like(pot),
            create_graph = True)[0]
        
        return drift
        
    def forward(self, x):
        return self.net(x)