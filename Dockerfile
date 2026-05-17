FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

WORKDIR /app

COPY serve.py ./
RUN chmod +x serve.py

# Pre-resolve script dependencies so the first run doesn't pay for it.
RUN uv run --script serve.py --help >/dev/null 2>&1 || true

COPY . .

CMD ["uv", "run", "--script", "serve.py", "9000"]
