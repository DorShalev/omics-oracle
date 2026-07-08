"""LLM wrapper with pluggable backends.

Backends, in auto-selection order:
  1. Claude CLI  — the local `claude` binary in print mode (subscription, no API key).
  2. LiteLLM     — OpenRouter/OpenAI etc. via an API key.
  3. none        — no model; every caller falls back to deterministic heuristics.

Auto rule: if ``OPENROUTER_API_KEY``/``OPENAI_API_KEY`` is set the user opted into the API,
so LiteLLM wins; otherwise the Claude CLI is used when the binary is present. Override with
``OMICS_ORACLE_LLM`` = ``auto`` | ``claude-cli`` | ``litellm`` | ``none``.

The model is used only to (a) turn a prompt into a structured spec and (b) phrase narrative
— never to produce numbers.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from typing import Optional

DEFAULT_LITELLM_MODEL = os.environ.get("OMICS_ORACLE_MODEL", "openrouter/anthropic/claude-3.5-sonnet")
# For the CLI, use its own short aliases ("sonnet"/"opus"/"haiku") or a full id; None => CLI default.
DEFAULT_CLI_MODEL = os.environ.get("OMICS_ORACLE_CLAUDE_MODEL")
CLI_TIMEOUT = float(os.environ.get("OMICS_ORACLE_CLI_TIMEOUT", "120"))


def _has_api_key() -> bool:
    return bool(os.environ.get("OPENROUTER_API_KEY") or os.environ.get("OPENAI_API_KEY"))


def _claude_binary() -> Optional[str]:
    return shutil.which("claude") or (
        os.path.expanduser("~/.local/bin/claude")
        if os.path.exists(os.path.expanduser("~/.local/bin/claude")) else None
    )


class LLM:
    """Public interface used across the engine: ``available``, ``complete``, ``complete_json``."""

    def __init__(self, model: Optional[str] = None, backend: Optional[str] = None):
        self.backend = (backend or os.environ.get("OMICS_ORACLE_LLM", "auto")).lower()
        self.litellm_model = model or DEFAULT_LITELLM_MODEL
        self.cli_model = DEFAULT_CLI_MODEL
        self._resolved = self._resolve_backend()

    def _resolve_backend(self) -> str:
        if self.backend in ("claude-cli", "cli"):
            return "claude-cli" if _claude_binary() else "none"
        if self.backend in ("litellm", "openrouter", "openai", "api"):
            return "litellm" if (_has_api_key() and _litellm_importable()) else "none"
        if self.backend == "none":
            return "none"
        # auto
        if _has_api_key() and _litellm_importable():
            return "litellm"
        if _claude_binary():
            return "claude-cli"
        return "none"

    @property
    def available(self) -> bool:
        return self._resolved != "none"

    @property
    def backend_name(self) -> str:
        return self._resolved

    def complete(self, system: str, user: str, temperature: float = 0.0,
                 max_tokens: int = 1024) -> Optional[str]:
        """Return raw text, or None if the model path is unavailable/fails."""
        if self._resolved == "claude-cli":
            return self._complete_cli(system, user)
        if self._resolved == "litellm":
            return self._complete_litellm(system, user, temperature, max_tokens)
        return None

    def complete_json(self, system: str, user: str) -> Optional[dict]:
        text = self.complete(system, user)
        return _extract_json(text) if text else None

    # -- backends -------------------------------------------------------------------------

    def _complete_cli(self, system: str, user: str) -> Optional[str]:
        binary = _claude_binary()
        if not binary:
            return None
        cmd = [binary, "-p", "--system-prompt", system]
        if self.cli_model:
            cmd += ["--model", self.cli_model]
        try:
            proc = subprocess.run(
                cmd, input=user, capture_output=True, text=True, timeout=CLI_TIMEOUT,
                # neutral cwd so the CLI doesn't scan this repo; it only needs to generate text
                cwd=os.path.expanduser("~"),
            )
        except (subprocess.TimeoutExpired, OSError):
            return None
        if proc.returncode != 0:
            return None
        out = (proc.stdout or "").strip()
        return out or None

    def _complete_litellm(self, system: str, user: str, temperature: float,
                          max_tokens: int) -> Optional[str]:
        try:
            import litellm
            resp = litellm.completion(
                model=self.litellm_model,
                messages=[{"role": "system", "content": system},
                          {"role": "user", "content": user}],
                temperature=temperature, max_tokens=max_tokens,
            )
            return resp["choices"][0]["message"]["content"]
        except Exception:
            return None


def _litellm_importable() -> bool:
    try:
        import litellm  # noqa: F401
        return True
    except Exception:
        return False


def _extract_json(text: str) -> Optional[dict]:
    # Prefer a fenced ```json block, else the first {...} span.
    fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.DOTALL)
    candidate = fence.group(1) if fence else None
    if candidate is None:
        brace = re.search(r"\{.*\}", text, re.DOTALL)
        candidate = brace.group(0) if brace else None
    if candidate is None:
        return None
    try:
        obj = json.loads(candidate)
        return obj if isinstance(obj, dict) else None
    except Exception:
        return None
