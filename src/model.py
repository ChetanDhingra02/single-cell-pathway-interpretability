import copy

import numpy as np
import torch
import torch.nn as nn


HIDDEN_UNITS = 256
LEARNING_RATE = 1e-3
BATCH_SIZE = 256
MAX_EPOCHS = 100
PATIENCE = 8
MIN_DELTA = 1e-4


class PathwayAutoencoder(nn.Module):

    def __init__(self, mask, nonnegative, hidden=HIDDEN_UNITS):
        super().__init__()

        n_latent, n_genes = mask.shape
        self.nonnegative = nonnegative

        self.encoder = nn.Sequential(
            nn.Linear(n_genes, hidden),
            nn.ReLU(),
            nn.Linear(hidden, n_latent)
        )

        if nonnegative:
            init = torch.rand(n_latent, n_genes) * 0.02
        else:
            init = torch.randn(n_latent, n_genes) * 0.01

        self.decoder_weight = nn.Parameter(init)
        self.decoder_bias = nn.Parameter(torch.zeros(n_genes))

        self.register_buffer(
            "mask",
            torch.as_tensor(mask, dtype=torch.float32)
        )


    @torch.no_grad()
    def clamp_nonnegative_(self):
        if self.nonnegative:
            self.decoder_weight.clamp_(min=0.0)


    def effective_weights(self):
        return self.decoder_weight * self.mask


    def encode(self, x):
        return self.encoder(x)


    def decode(self, z):
        return z @ self.effective_weights() + self.decoder_bias


    def forward(self, x):
        z = self.encode(x)
        return self.decode(z), z


def train_autoencoder(
    X_train,
    X_val,
    mask,
    nonnegative,
    seed,
    max_epochs=MAX_EPOCHS
):
    torch.manual_seed(seed)

    model = PathwayAutoencoder(
        mask,
        nonnegative
    )

    optimiser = torch.optim.Adam(
        model.parameters(),
        lr=LEARNING_RATE
    )

    loss_fn = nn.MSELoss()

    train = torch.from_numpy(X_train)
    val = torch.from_numpy(X_val)

    generator = torch.Generator().manual_seed(seed)

    history = []
    best_loss = np.inf
    best_state = None
    best_epoch = 0
    waited = 0

    for epoch in range(1, max_epochs + 1):

        model.train()

        order = torch.randperm(
            len(train),
            generator=generator
        )

        for start in range(0, len(train), BATCH_SIZE):

            batch = train[
                order[start:start + BATCH_SIZE]
            ]

            loss = loss_fn(
                model(batch)[0],
                batch
            )

            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            model.clamp_nonnegative_()

        model.eval()

        with torch.no_grad():
            train_loss = loss_fn(
                model(train)[0],
                train
            ).item()

            val_loss = loss_fn(
                model(val)[0],
                val
            ).item()

        assert np.isfinite(train_loss)
        assert np.isfinite(val_loss)

        history.append({
            "epoch": epoch,
            "train_loss": train_loss,
            "val_loss": val_loss
        })

        if val_loss < best_loss - MIN_DELTA:

            best_loss = val_loss
            best_epoch = epoch
            waited = 0
            best_state = copy.deepcopy(
                model.state_dict()
            )

        else:

            waited += 1

            if waited >= PATIENCE:
                break

    model.load_state_dict(best_state)
    model.eval()

    return model, history, best_epoch, best_loss


def train_fixed_epochs(X_train, mask, nonnegative, seed, n_epochs):
    # Final training uses the epoch count chosen from donor CV
    torch.manual_seed(seed)

    model = PathwayAutoencoder(mask, nonnegative)
    optimiser = torch.optim.Adam(
        model.parameters(),
        lr=LEARNING_RATE
    )
    loss_fn = nn.MSELoss()

    train = torch.from_numpy(X_train)
    generator = torch.Generator().manual_seed(seed)

    history = []

    for epoch in range(1, n_epochs + 1):

        model.train()

        order = torch.randperm(
            len(train),
            generator=generator
        )

        for start in range(0, len(train), BATCH_SIZE):

            batch = train[
                order[start:start + BATCH_SIZE]
            ]

            loss = loss_fn(
                model(batch)[0],
                batch
            )

            optimiser.zero_grad()
            loss.backward()
            optimiser.step()

            model.clamp_nonnegative_()

        model.eval()

        with torch.no_grad():
            train_loss = loss_fn(
                model(train)[0],
                train
            ).item()

        history.append({
            "epoch": epoch,
            "train_loss": train_loss
        })

    return model, history


def latent_scores(model, X, batch=4096):

    scores = []

    with torch.no_grad():

        for start in range(0, len(X), batch):

            scores.append(
                model.encode(
                    torch.from_numpy(
                        X[start:start + batch]
                    )
                ).numpy()
            )

    return np.vstack(scores)


def reconstruction_mse(model, X):

    with torch.no_grad():

        x = torch.from_numpy(X)

        return float(
            torch.mean(
                (model(x)[0] - x) ** 2
            )
        )