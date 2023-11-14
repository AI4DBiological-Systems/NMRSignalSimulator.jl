using BenchmarkTools
using Test
using SparseArrays
using DataDeps, Tar, CodecZlib

using LinearAlgebra
using FFTW

import Random

import PythonPlot as PLT
import PublicationDatasets as DS
import NMRDataSetup as DSU

using BenchmarkTools
import Statistics

using Revise

import NMRSignalSimulator
SIG = NMRSignalSimulator

NMRHamiltonian = NMRSignalSimulator.NMRHamiltonian
JSON3 = NMRHamiltonian.JSON3

HAM = NMRHamiltonian