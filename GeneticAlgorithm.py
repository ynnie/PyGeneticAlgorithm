import random
import numpy
import math

class GeneticAlgorithm:
# Genetic Algorithm
    
    # Constructer
    def __init__(self, maxGen, popGeneration, pCross, pMutation, nGenes, parametersRange, fitnessFunc):
        self.maxGenNum       = maxGen           # Maxium generation number
        self.popGeneration   = popGeneration    # Population of one generation
        self.pCross          = pCross           # Probability of cross interchange
        self.pMutation 		 = pMutation        # Probability of mutation
        self.nGenes    		 = nGenes           # Number of genes(parameters)
        self.parametersRange = parametersRange  # A matrix of parameters' ranges, [lower, upper, accuracy]
        self.fitnessFunc     = fitnessFunc      # Fitness function, parameters as input, fitness as return

        self.generation      = []
        self.fitnessList     = []
        self.bestFitness     = 0.0
        self.dnaIdx          = []
        self.geneLength      = []
        idx = 0
        for i in range(nGenes):
            lowerLimit  = parametersRange[i][0]
            upperLimit  = parametersRange[i][1]
            accuracy    = parametersRange[i][2]
            parRange    = upperLimit - lowerLimit
            geneLength  = math.ceil(math.log2(parRange/accuracy))+1
            self.geneLength.append(geneLength)
            startIndex  = idx
            idx         = idx + geneLength
            endIndex    = idx - 1
            self.dnaIdx.append([startIndex, endIndex])
        self.dnaLength  = idx
        self.bestIndividual = list('-'*self.dnaLength)
        return

    # Display information about current object
    def info(self):
        print('Genetic Algorithm:')
        print('Max generation number:', self.maxGenNum)
        print('Population of one generation:', self.popGeneration)
        print('Probability of mutation:', self.pMutation)
        print('Probability of cross:', self.pCross)
        print('Number of genes(parameters):', self.nGenes)
        print('Range of parameters: [lower, upper, accuracy] ', self.parametersRange)
        return
          
    # Encode parameters into DNA    
    def encode(self, parameters):
        dna = ''
        for i in range(self.nGenes):
            geneSequenceDigit = (parameters[i]-self.parametersRange[i][0])/(self.parametersRange[i][1]-self.parametersRange[i][0])*(2**self.geneLength[i]-1)
            geneSequenceDigit = int(round(geneSequenceDigit))
            geneSequenceBin   = self.int2bin(geneSequenceDigit, self.geneLength[i])
            dna = dna + geneSequenceBin 
        dna = list(dna)  # Trun string to list
        return dna
    
    # Decode DNA to parameters
    def decode(self, dna):
        dna = ''.join(dna)  # Trun list to string
        parameters = []
        for i in range(self.nGenes):
            geneSequenceBin   = dna[self.dnaIdx[i][0]:self.dnaIdx[i][1]+1]
            geneSequenceDigit = self.bin2int(geneSequenceBin)
            parameterI = geneSequenceDigit/(2**self.geneLength[i]-1)*(self.parametersRange[i][1]-self.parametersRange[i][0])+self.parametersRange[i][0]
            parameters.append(parameterI)
        return parameters
    
    # Create random DNA for initially generation
    def randomDNA(self):
        dna = ''
        for i in range(self.nGenes):
            geneSequenceDigit = random.randint(0, 2**self.geneLength[i])
            geneSequenceBin   = self.int2bin(geneSequenceDigit, self.geneLength[i])
            dna = dna + geneSequenceBin
        dna = list(dna) # Trun string to list
        return dna
    
    # Generate initially generation
    def iniGeneration(self):
        for i in range(self.popGeneration):
            self.generation.append(self.randomDNA())
        self.fitness()
        self.bestFitness = max(self.fitnessList)
        idx = self.fitnessList.index(self.bestFitness)
        self.bestIndividual = self.generation[idx].copy()
        return
    
    # Mutation
    def mutation(self):
        for i in range(self.popGeneration):
            ret = random.random()
            while ret == 0:
                ret = random.random()
            if ret < self.pMutation:
                mutationPoint = random.randint(0,self.dnaLength-1)
                if self.generation[i][mutationPoint] == '0':
                    self.generation[i][mutationPoint] = '1'    
                else:
                    self.generation[i][mutationPoint] = '0'
        return
    
    # Cross interchange
    def cross(self):
        for i in range(self.popGeneration):
            ret = random.random()
            while ret == 0:
                ret = random.random()
            if ret < self.pCross:
                pair = random.randint(0,self.popGeneration-1)
                while pair == i:
                    pair = random.randint(0,self.popGeneration-1)
                crossPiont = random.randint(0,self.dnaLength-1)
                self.generation[i][crossPiont:-1] = self.generation[pair][crossPiont:-1]
        return
    
    # Select individuals
    def select(self):
        nCandidates = 2
        generation  = []
        for i in range(self.popGeneration):
            candidates = random.choices(self.generation, k=nCandidates)
            fitnessList = self.fitness(candidates=candidates)
            bestFitness = max(fitnessList)
            idx = fitnessList.index(bestFitness)
            generation.append(candidates[idx].copy())
        self.generation = generation
        return
    
    # Elitist preservation
    def elitistPreservation(self):
        bestFitness = max(self.fitnessList)
        if bestFitness > self.bestFitness:
            idx = self.fitnessList.index(bestFitness)
            self.bestFitness    = bestFitness
            self.bestIndividual = self.generation[idx].copy()
        else:
            worstIndividual = min(self.fitnessList)
            idx = self.fitnessList.index(worstIndividual)
            self.generation[idx]  = self.bestIndividual.copy()
            self.fitnessList[idx] = self.bestFitness
        return
    
    # Fitness method for selecting individuals
    def fitness(self, candidates='None'):
        if candidates=='None':
            self.fitnessList = []
            for dna in self.generation:
                parameters = self.decode(dna)
                self.fitnessList.append(self.fitnessFunc(parameters))
            return
        else:
            fitnessList = []
            for dna in candidates:
                parameters = self.decode(dna)
                fitnessList.append(self.fitnessFunc(parameters))
            return fitnessList

    
    # Iterate 'maxGen' generations to find the best individual
    def run(self):
        self.iniGeneration()
        print('Genetic Algorithm is running:')
        for i in range(1, self.maxGenNum+1):
            self.select()
            self.cross()
            self.mutation()
            self.fitness()
            self.elitistPreservation()
            
            if i%1 == 0:
                meanFitness = numpy.mean(self.fitnessList)
                print('N:', i,'Best:', ''.join(self.bestIndividual), '%.3f' %self.bestFitness, 'Mean:%.3f' %meanFitness)
        return

    def result(self):
        meanFitness = numpy.mean(self.fitnessList)
        print('Best:', ''.join(self.bestIndividual), '%.3f' %self.bestFitness, 'Mean:%.3f' %meanFitness)

    # Returns the binary string of integer n, using count number of digits
    def int2bin(self, n, count=32):
        binStr = "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])
        return binStr

    # Returns digit integer of binary string
    def bin2int(self, n):
        ret = int(n,2)
        return ret