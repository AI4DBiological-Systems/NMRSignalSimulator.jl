using BenchmarkTools
using Test
using SparseArrays
using DataDeps, Tar, CodecZlib

using LinearAlgebra

import Random

import PythonPlot

using BenchmarkTools
import Statistics

using Revise

import NMRSignalSimulator
SIG = NMRSignalSimulator

NMRHamiltonian = NMRSignalSimulator.NMRHamiltonian
JSON3 = NMRHamiltonian.JSON3

HAM = NMRHamiltonian