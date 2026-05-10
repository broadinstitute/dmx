---
name: getting-started
description: >-
  Walk a first-time dmx user from a fresh clone to a running marimo kernel
  with agent composition enabled. Trigger when the user says "help me get
  started", "onboard me", "set me up", "I'm new to dmx", or asks a
  DepMap/Breadbox composition question before a marimo kernel is running.
---

# Getting started with dmx

Your job: get this user from a cold clone to a live marimo kernel, then hand
off to the `compose-notebook` skill for the actual composition.

## Setup flow

### 1. Verify uv is installed

Run `uv --version`. If it fails, tell the user to run:

    curl -LsSf https://astral.sh/uv/install.sh | sh

Then have them source their shell profile or open a new terminal. Re-check
`uv --version`.

### 2. Install marimo skills

Install the upstream marimo authoring and live-kernel skills globally:

    npx skills add marimo-team/skills -g --agent codex -y
    npx skills add marimo-team/marimo-pair -g --agent codex -y

If the session cannot see marimo-pair after install, ask the user to restart
Codex and run `getting-started` again.

### 3. Launch the orientation notebook

From the dmx repo root:

    PORT=$(python -c "import socket; s=socket.socket(); s.bind(('127.0.0.1',0)); print(s.getsockname()[1])")
    env -u PYTHONPATH uvx marimo edit --sandbox --headless --no-token \
        --port $PORT notebooks/nb01_orientation.py

Run it in the background. Verify with:

    curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:$PORT/

Expect HTTP 200. Tell the user the URL.

### 4. Hand off

Once the kernel is live, ask what they want to explore and use
`compose-notebook`. The typical first move is dataset discovery (`nb02`) or
a gene dependency profile (`nb03`).

## Gotchas

- Use `--sandbox`; the notebooks rely on PEP 723 headers for dependencies.
- Use `env -u PYTHONPATH` on machines where Nix or other shells inject paths.
- Breadbox public reads require no token, but the service is remote. Network
  or API failures are runtime failures, not notebook syntax failures.
