# Load Packages
using Ipopt, JuMP, DataFrames

# Read in Data
data = readtable("Data.csv")
IntData = readtable("IntData.csv")

# Creat JuMP Model
M = Model(solver = IpoptSolver())


# PARAMETERS #

## Production block

QCapital0 = data[:QCapital0]              # inital capital demand
CapSupply0 = sum(QCapital0)               # capital endowment
QLabor0  = data[:QLabor0]                 # initial labor demand
LaborSupply0 = sum(QLabor0)               # labor endowment
PCapital0 = 1.                            # initial capital price
PLabor0 = 1.                              # initial labor price (wages)
PComProd0  = data[:PComProd0]             # initial commodity price (price of domestically-produced commodities)
IntInputs0 = [IntData[:AgCom] IntData[:IndCom] IntData[:ServCom]]'
                                          # intermediate inputs
QComProd0 = (PCapital0.*QCapital0)./PComProd0 + (PLabor0.*QLabor0)./PComProd0 +        [sum(IntInputs0[:,1]),sum(IntInputs0[:,2]),sum(IntInputs0[:,3])].*PComProd0
                                          # initial commodity production level
intcoeff = IntInputs0./QComProd0'         # technical coefficients
elasF  = data[:elasF]                     # elasticity of subsistition between factors in production function
distF = 1./(1.+((PLabor0/PCapital0)*((QCapital0./QLabor0).^(-1./elasF))))
                                          # distribution parameter
eF = QComProd0 ./ (distF.*QCapital0.^((elasF-1)./elasF)+(1-distF).*QLabor0.^((elasF-1)./elasF)).^(elasF./(elasF-1))
                                          # efficiency parameter in production function

## Consumption block

QComCons0  = data[:QComCons0]             # household demand for commodities
Income0 = PCapital0*CapSupply0 + PLabor0*LaborSupply0
                                          # household's total income
frisch = -1.1			                        # Frisch parameter
elasI  = data[:elasI]                     # income elasticities of demand for commodities
sharesULES = (elasI.*((PComProd0.*QComCons0)./Income0))./sum(elasI.*((PComProd0.*QComCons0)./Income0))
                                          # marginal budget shares of utility function
subsist = QComCons0 + (sharesULES.*Income0)./(PComProd0.*frisch)
                                          # subsistence level
Utility0 = prod((QComCons0 - subsist).^sharesULES)
                                          # utility level of household


# VARIABLES #

## Define indices

sector = [1,2,3] # Dimension of sectors
com = [1,2,3]    # Dimension of commodities

## Production block

@variables M begin
  QCapital[i = sector], (start = QCapital0[i])
  PCapital, (start = PCapital0)
  QLabor[i = sector], (start = QLabor0[i])
  PLabor == PLabor0                               # numeraire
  PComProd[i = sector], (start = PComProd0[i])
  CapSupply == CapSupply0                         # fixed capital supply
  LaborSupply == LaborSupply0                     # fixed labor supply
  QComProd[i = sector], (start = QComProd0[i])
end

## Consumption block

@variables M begin
  QComCons[i = sector], (start = QComCons0[i])
  Income, (start = Income0)
  Utility, (start = Utility0)
end


# EQUATIONS/CONSTRAINTS #

## Production block

@NLconstraints M begin
  EQ_QCapital[i = sector], QCapital[i] == distF[i]^elasF[i] * PCapital^(-elasF[i]) * (distF[i]^elasF[i] * PCapital^(1-elasF[i]) + (1-distF[i])^elasF[i] * PLabor^(1-elasF[i]))^(elasF[i]/(1-elasF[i])) * (QComProd[i]/eF[i])
                                                                       # firm demand for capital
  EQ_QLabor[i = sector], QLabor[i] == (1-distF[i])^elasF[i] * PLabor^(-elasF[i]) * (distF[i]^elasF[i] * PCapital^(1-elasF[i]) + (1-distF[i])^elasF[i] * PLabor^(1-elasF[i]))^(elasF[i]/(1-elasF[i])) * (QComProd[i]/eF[i])
                                                                       # firm demand for labor
  EQ_ZPC[i = sector], PComProd[i] * QComProd[i] == PCapital*QCapital[i] + PLabor*QLabor[i] + sum(intcoeff[l,i] * PComProd[l] for l in sector)*QComProd[i]
                                                                       # zero-profit condition
  EQ_CapSupply, sum(QCapital[j] for j in sector) == CapSupply           # capital market clearing
  EQ_LaborSupply, sum(QLabor[j] for j in sector) == LaborSupply         # labor market clearing
end

## Consumption block

@NLconstraints M begin
  EQ_QComCons[i = sector], QComCons[i]  == subsist[i]  + (sharesULES[i] /PComProd[i]) * (Income - sum(PComProd[j]  * subsist[j] for j in sector ))
                                                                       # consumer consumption
  EQ_Income, Income == PCapital*CapSupply + PLabor*LaborSupply	         # income balance
  EQUtility, Utility == prod((QComCons[j]  - subsist[j])^sharesULES[j] for j in sector)
                                                                       # household utility
  EQ_QComProd[i = sector], QComProd[i]  == sum(intcoeff[i,j]  * QComProd[j] for j in sector) + QComCons[i]	                                                         # market clearing consumption
end


# Model Solver

@NLobjective(M, Max, 0)
@time solve(M)
