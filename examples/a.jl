using BenchmarkTools

using LinearAlgebra
using FFTW

import Random

import PythonPlot as PLT
import PublicationDatasets as DS
import NMRDataSetup as DSU

import Statistics

using Revise

import NMRSignalSimulator
SIG = NMRSignalSimulator

NMRHamiltonian = NMRSignalSimulator.NMRHamiltonian
JSON3 = NMRHamiltonian.JSON3

HAM = NMRHamiltonian