"""
Unit and regression test for the photonic_circuit_solver package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import photonic_circuit_solver


def test_photonic_circuit_solver_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "photonic_circuit_solver" in sys.modules
