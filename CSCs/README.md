# In silico mesoscopic model to simulate tumor growth and vasculature based on cancer stem cells hypothesis

This is a work based on the mesoscopic model of Jiménez-Sánchez et al., (2021). We extended it to the Cancer Stem Cells (CSCs) hypothesis. Here, the details are not discussed, but we rather give a quick review of the packages needed to run the code correctly.

## 1. Julia Installation and packages installation

Please go to the READme of the main repository, a very detailed guideline written by Jiménez-Sánchez et al. is found there.
But here, we would like to add a small detail for installing packages when running Julia in Windows Powershell. This is the interface you will find:

```
julia>
```
Here just type `]`, this following interface shows up:

```
(@v1.12) pkg>
```
The parenthesis' content depends on your `Julia` version. Here, you can add or see the packages available in your "environment". To add a package:

```
(@v1.12) pkg> add Optuna
```

To see the packages available:

```
(@v1.12) pkg> status
(@v1.12) pkg> status
Status `your_path2Julia\.julia\environments\v1.12\Project.toml`
  [31c24e10] Distributions v0.25.123
  [a5d0552b] Optuna v0.2.1
  [6099a3de] PythonCall v0.9.31
⌃ [295af30f] Revise v3.14.0
Info Packages marked with ⌃ have new versions available and may be upgradable.
```

## 2. Julia packages requirements
Basically all packages requiered to run the main.jl are installed after Julia's instalation. Only an exclusive package is needed to run `parameter_optimization.jl`: PythonCall. You can install it easily:

```
(@v1.12) pkg> add PythonCall
```

This is called in Script as:

```
26 using PythonCall
```

This enables the usage of some Python Modules, which is very convinent and useful when dealing with machine-learning problems. Here, in the Script `parameter_optimization.jl`, `Optuna` module is called from Python:

```
36 pyimport ("Optuna")
```
Since its usage is much simpler than the `Optuna` package available in Julia. However, there are some inconvinience as some syntax have to be changed from Julia's one. For example:
```
280 CSC_rate = pyconvert(Float64, trial.params["CSC_rate"])
```

## 3. Python Packages requirements
Python is only used for the parameter optimization's results analysis. The most basic ones for machine-learning and image plotting were used, and all of them can be easily installed via `pip` or `conda`:

```
pip install numpy
pip install pandas
pip install matplotlib
pip install scipy
```

## 4. Contacts and the objective of this project.
This project was developed mainly due to the inscription of `Math4MP`, an extra activity developed by UPM, with the main objective of using Mathematics in different disciplines. But also, this project used for competing in  `EELISA Scientific Students Competition 2026`. 

Contacts:

Patricia Martín Reguera. Email: patricia.martinreg@alumnos.upm.es // Github:https://github.com/patmareg

Yixi Zhang. Email: yixi.z@alumnos.upm.es // Github: https://github.com/YixiZ04
 



