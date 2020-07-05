# PyGeneticAlgorithm
A sample Genetic Algorithm implementation in Python.

## How to use:

To use this project, you need to clone this module to your project folder.

1. Import the GeneticAlgorithm module by `import GeneticAlgorithm as GA`;

2. Set your own fitness function, the fitness function must accept the parameters that you want to optimize and return a fitness score;

3. Create a GeneticAlgorithm subject : `ga = GA.GeneticAlgorithm(maxGen, popGeneration, pCross, pMutation, nGenes, parametersRange, fitness)`;

   ```python
   maxGen # How many generation
   popGeneration # How many individuals in one geration
   pCross # The possibility of cross
   pMutation # The possibility of mutation
   nGenes  # How many genes (parameters to be optimized)
   parametersRange # The ranges of parameters, e.g. [[-1, 2, 0.01]，[-2, 3, 0.1]]
   fitness # The fitness function you defined
   ```

4. Run: `ga.run()`

5. Get result: `ga.result()`

## Example: 

```python
import GeneticAlgorithm as GA # Import GeneticAlgorithm module
import math

# Set fitness function, par is the parameter list，ret is the fitness score
def fitness(par):
    x = par[0]
    ret = x*math.sin(10*math.pi*x)+2
    return ret

maxGen = 500 # How many generation
popGeneration = 50 # How many individuals in one geration
pCross = 0.1 # The possibility of cross
pMutation = 0.05 # The possibility of mutation
nGenes = 1 # How many genes (parameters to be optimized)
parametersRange = [[-1, 2, 0.01]] # The ranges of parameters, e.g. [[-1, 2, 0.01]，[-2, 3, 0.1]]

# Create a GA subject
ga = GA.GeneticAlgorithm(maxGen, popGeneration, pCross, pMutation, nGenes, parametersRange, fitness)
# Display subject information
ga.info()
# Run Genetic Algorithm
ga.run()
# Show result
ga.result()
```