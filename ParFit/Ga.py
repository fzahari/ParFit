#!/usr/bin/env python

import random,numpy
from deap import algorithms, base, creator, tools

def varAnd(population, toolbox, cxpb, mutpb):

    offspring = [toolbox.clone(ind) for ind in population]

    # Apply crossover and mutation on the offspring
    for i in range(1, len(offspring), 2):
        if random.random() < cxpb:
            offspring[i - 1], offspring[i] = toolbox.mate(offspring[i - 1], offspring[i])
            del offspring[i - 1].fitness.values, offspring[i].fitness.values

    for i in range(len(offspring)):
        if random.random() < mutpb:
            offspring[i], = toolbox.mutate(offspring[i])
            del offspring[i].fitness.values

    return offspring

def eaSimple(population, toolbox, cxpb, mutpb, ngen, stats=None,
             halloffame=None, verbose=__debug__):

    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print logbook.stream

    flag=0
    best_ind=halloffame[0] 
 
    # Begin the generational process
    for gen in range(1, ngen + 1):
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population))

        # Vary the pool of individuals
        offspring = varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        if numpy.allclose(numpy.array(best_ind),numpy.array(halloffame[0])):
           flag+=1 
        else:
           best_ind=halloffame[0] 
           flag=0

        # Replace the current population by the offspring
        population[:] = offspring

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print logbook.stream

        #print flag, best_ind, halloffame[0]
        if flag==4: break

    return population, logbook

def run_ga(engine_rmse,np,ngen):

   creator.create("FitnessMax",base.Fitness,weights=(-1.0,))
   creator.create("Individual",numpy.ndarray,fitness=creator.FitnessMax)

   toolbox=base.Toolbox()
   toolbox.register("attr_float",numpy.random.random)
   toolbox.register("individual",tools.initRepeat,creator.Individual,toolbox.attr_float,n=np)
   toolbox.register("population",tools.initRepeat,list,toolbox.individual)
   toolbox.register("evaluate",engine_rmse)
 
   def cxComb(ind1,ind2):
      if numpy.random.randint(2)==0:
         tools.cxOnePoint(ind1,ind2)
      else:
         tools.cxBlend(ind1,ind2,0.0)
      return ind1,ind2

   toolbox.register("mate",cxComb)
   toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=1.)

   def checkBounds(min, max):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in xrange(len(child)):
                    if child[i] > max:
                        child[i] = max
                    elif child[i] < min:
                        child[i] = min
            return offspring
        return wrapper
    return decorator

   MIN=-3.0
   MAX=3.0

   toolbox.decorate("mate", checkBounds(MIN, MAX))
   toolbox.decorate("mutate", checkBounds(MIN, MAX))

   toolbox.register("select",tools.selTournament,tournsize=3)

   stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
   stats.register("avg", numpy.mean)
   stats.register("std", numpy.std)
   stats.register("min", numpy.min)
   stats.register("max", numpy.max)

   hof = tools.HallOfFame(np, similar=numpy.allclose)

   pop=toolbox.population(n=50)
   eaSimple(pop,toolbox,cxpb=0.3,mutpb=0.05,ngen=ngen,stats=stats,halloffame=hof,verbose=True)
   print "The best rmse is", hof[0].fitness.values[0],"corresponding to the parameters",numpy.array(hof[0])

   #return hof[0:int(np/20)]
   return hof
