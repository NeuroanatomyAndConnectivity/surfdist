"""Smoke test: package and its declared submodules import cleanly."""

import importlib

import surfdist


def test_submodules_importable():
    for name in surfdist.__all__:
        mod = importlib.import_module(f'surfdist.{name}')
        assert mod is not None, f'surfdist.{name} failed to import'
