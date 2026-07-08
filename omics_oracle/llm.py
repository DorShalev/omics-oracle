"""Thin LLM wrapper over LiteLLM/OpenRouter.

The engine must run with *no* key: every call site has a deterministic fallback, and
`LLM.available` tells callers whether the model path is live. The model is used only for
(a) turning a prompt into a structured spec and (b) phrasing narrative — never to produce
numbers.
"""
from __future__ import annotations

import json
import os
import re
from typing import Optional


DEFAULT_MODEL = os.environ.get("OMICS_ORACLE_MODEL", "openrouter/anthropic/claude-3.5-sonnet")


class LLM:
    def __init__(self, model: Optional[str] = None):
        self.model = model or DEFAULT_MODEL
        # OpenRouter key is the primary path; LiteLLM also honors provider-native keys.
        self._key = os.environ.get("OPENROUTER_API_KEY") or os.environ.get("OPENAI_API_KEY")

    @property
    def available(self) -> bool:
        if not self._key:
            return False
        try:
            import litellm  # noqa: F401
        except Exception:
            return False
        return True

    def complete(self, system: str, user: str, temperature: float = 0.0,
                 max_tokens: int = 1024) -> Optional[str]:
        """Return raw text, or None if the model path is unavailable/fails."""
        if not self.available:
            return None
        try:
            import litellm
            resp = litellm.completion(
                model=self.model,
                messages=[
                    {"role": "system", "content": system},
                    {"role": "user", "content": user},
                ],
                temperature=temperature,
                max_tokens=max_tokens,
            )
            return resp["choices"][0]["message"]["content"]
        except Exception:
            return None

    def complete_json(self, system: str, user: str) -> Optional[dict]:
        """Complete and parse a JSON object out of the response, tolerating code fences."""
        text = self.complete(system, user)
        if not text:
            return None
        return _extract_json(text)


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
