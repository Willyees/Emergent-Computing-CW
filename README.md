# Coursework 2 Emergent computing optimization - Napier University 2019/2020
"Team Pursuit Track Cycling Simulator"

## Aim
A program is provided in which a cycling race in a velodrome is simulated. The aim of the cw is to improve the time taken to complete the race by applying Evolutionary Algorithms and find out the best strategy in transition between first cyclist and the pace to be utilized during each race sections.

![velodrome](https://github.com/Willyees/Emergent-Computing-CW/blob/assets/assets/velodrome.jpg)

## Methodology

![initial_times](https://github.com/Willyees/Emergent-Computing-CW/blob/assets/assets/initial_times.png)

Multiple GitHub branches have been created in which different Evolutionary Algorithms have been attempted to solve the problem:
hill climber, island, restricted mating, simulated annealing and speciation.
Additionally, it has been explored the effect of modifying default parameters, understanding their effect on the race time and on the algorithm learning-speed. 
Some variants of main components of the EAs have been experimented (for example the use of a tournament selection rather than a best children chosen)
The fitness function has been explored, aiming to obtain a good optimum in a mid-long time algorithm run.


## Observations
By observing the algorithm running it was noticed a problem of early convergence, leading to one of the biggest problem of the algorithm: local optima. This was noticed to be a product of a very static population with many similar (even same) chromosomes. Variations of inner functionalities have been explored to find out the best to be applied.
Additionally, it has been quite challenging to improve the fitness function, because it was found that the search space is very vast but the number of correct chromosomes that would lead to a finishing race is very slim compared to the total. Without a fitness function that would direct the search towards correct chromosomes, it would take a very long time to obtain correct solutions just by mutating (very unlikely since search space is vast).
After implementing these algorithms it was understood that the most suitable are the ones that can easily work with a rough search space. The hill climber is not very suitable and the simulated annealing works slightly better. 
The EA and the Island are the most suitable because there were implemented functions to diversify population and escape optimums.


More in detail info, observations and procedures followed in improving the race time can be found in the [report](./report.pdf)

![final_times](https://github.com/Willyees/Emergent-Computing-CW/blob/assets/assets/final_times.png)

======================================

Creator of the simulator: Markus Wagner (wagner@acrocon.com)
Used in: http://cs.adelaide.edu.au/~optlog/research/cycling.php
Paper:   Evolving Pacing Strategies for Team Pursuit Track Cycling
         Markus Wagner, Jareth Day, Diora Jordan, Trent Kroeger, and Frank Neumann
         arxiv:    http://arxiv.org/abs/1104.0775
         Springer: http://link.springer.com/chapter/10.1007/978-1-4614-6322-1_4
         slides (MIC 2011): https://cs.adelaide.edu.au/users/honours/ec/2011-s2/cycling.pdf

Solution representation:
- Women's event (3 cyclists, all have to reach the finish): 23 integers for the pacing part, 22 booleans for the transition part
