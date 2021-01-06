traveling salesman problem (TSP) using Genetic-Algorithms

## Introduction

The traveling salesman problem (TSP) is a typical example of a very hard combinatorial optimization problem. The problem is to find the shortest tour that passes through each vertex in a given graph exactly once.  We will be use the GA approach to solve this problem. This is described in the following algorithm. 
Genetic Algorithm Psedocode:
Initialize population 
while evolution do 
Begin
 	choose parentl and parent2; { selection }
 	offspring := combination (parentl, parent2); 
{crossover} optimize-local (offspring); 
{mutation} ifsuited (offspring)  survival of the fittest 
then perform mutation through swap
 end;
In the following sections we will discuss each step of the above algorithm in detail. 
1. Initialization
This is done by declaring a 2D array and assigning each member a random value. This results in different combinations of tours Generated. . 
## 2. Evaluation and selection:
GAs use a selection mechanism to select individuals from the population to insert into a mating pool. Individuals from the mating pool are used to generate new offspring , with the resulting offspring forming the basis of the next genera tion. As the individuals in the mating pool are the ones whose genes are inherited by the next generation, it is desirable that the mating pool be comprised of "good" individuals. A selection mechanism in GAs is simply a process that favors the selection of better individuals in the populat ion for the mating pool.
To improve the solutions, each combination in the population is evaluated using a measure of fitness. Using the fitness function we can find the best Tour(combination) that corresponds to the shortest path. 
As TSP is a minimization problem so to convert it into maximization problem we have considered fitness function f(x) =1/d, where d calculates cost (or distance/length) of the tour represented by a chromosome. 

### Normalization:
To normalize the results we multiplied it a with a large constant(e.g 10^6). This results the fitness to mostly lie in a range of [0,30]. The fitness function that characterizes each route represents the total length of the route from the first to the last gene (city) moving according to the order of the genes in the route.
Distance Calculation:
 If the cities are represented with x and y coordinates in 2D coordinate system, then we calculate the distance between them according the equation: 
 
To avoid recalculating the distance each time we evaluate the fitness of a tour, we implemented a n x n table which store the distance from every single city to every other city. 
### 2.1 Fitness proportionate selection using Roulette-Wheel:
Roulette wheel selection is a frequently used selection operator in implementation of GA. In this method each individual is assigned a slice of a circular “roulette wheel”, the size of the slice being proportional to the individual’s fitness. The wheel is spun N times, where N is number of individuals in the population. On each spin, the individual under the wheel’s marker is selected to be in the pool of parents for the next generation.
 If   is the fitness of individual   in the population, its probability of being selected is given by the following equation where   is the number of individuals in the population.
 
### PseduoCode:
for all members of population
    sum += fitness of this individual
end for
for all members of population
    probability = sum of probabilities + (fitness / sum)
    sum of probabilities += probability
end for
This results in each member of the population being assigned a range of probability. After this step we can perform selection:
        number = Random between 0 and 1
        for all members of population
 if number > probability but less than next probability 
                then you have been selected
### 2.2 Tournament Selection:
Tournament selection is a method of selecting an individual from a population of individuals in a genetic algorithm. Tournament selection involves running several "tournaments" among a few individuals chosen at random from the population. The winner of each tournament (the one with the best fitness) is selected for crossover. Selection pressure is easily adjusted by changing the tournament size. If the tournament size is larger, weak individuals have a smaller chance to be selected. 

##3. Generation of offspring population:
To create the next generation, a new set of population called offspring is formed by the execution of genetic operators such as selection, Crossover and mutation. In particular, the crossover operator acts as a main operator and exercises a great influence on the performance of the GA approach, while the mutation operator acts as a background operator.
###3.1 Crossover:
Many crossover techniques exist for organisms which use different data structures to store themselves. , a direct swap may not be possible. One such case is when the chromosome is an ordered list, such as an ordered list of the cities to be travelled for the traveling salesman problem. There are many crossover methods for ordered chromosomes. 
Since Crossover is considered to a very important step in GA, We did an extensive literature review of all the popular and effective techniques that exist.
P. LARRANAGA, C.M.H, R.H. and Kylie Bryant have performed a good review and comparison of all the existing methods. Notable Methods Include: 
 

We experimented with the following two techniques for implementing the cross over function:
### 3.1.1. PMX crossover Method : 
Explained by Gokturk  Ucoluk, PMX crossover works as follows:
 Given two parents s and t, PMX randomly picks a crossover point – like 1-point crossover. The child is then constructed in the following way. Starting with a copy of s, the positions between the crossover points are, one by one, set to the values of t in these positions. To keep the string a valid chromosome the cities in these positions are not just overwritten. To set position p to city c, the city in position p and city c swap positions. Below you see an example of this coding and special crossover technique for two sample permutations: 5, 7, 1, 3, 6, 4, 2 and 4, 6, 2, 7, 3, 1, 5.
 
### 3.1.2 Moon Crossover:
Proposed by C. Moon et al. Moon crossover is novel method for crossover operation. The procedure of the moon crossover operator is described in Fig. 4. 
Example:  Let take an example to further understand this method. In Fig. 4, <osp, gi> is the concatenation operator to add gi after the substring osp. An example of the moon crossover is illustrated in Fig. 5. Suppose that two chromosomes are pa =[5 1 7 2 4 6] and pb is  [3 6 1 4 2 5 7]. First, select the substring from pa at random. In this example, the substring is selected as osp is [7 2]. Then, we can obtain sub pb [3 6 1 4 5] from pb. Next, g2 is 1 and q1 is 3 because i<-3 -1 and k <-0 + 1. The offspring becomes osp [1 7 2 3]. In the same way, add g1; q2, and the offspring becomes osp [5 1 7 2 3]. Now the next is g7 is 3 and q3 is  1 and the cities have already appeared in the offspring, so we cannot add these cities into the offspring osp. Finally, q4 is picked and the offspring becomes osp = [5 1 7 2 3 6 4]. Using the same procedure we can produce the second offspring as [3 6 1 4 5 7 2]. 

 
### 3.2 Mutation:
The swap mutation operates by selecting the two genes within a tour(population member) and then swaping their values. The probability of this operation happening is kept very low. 
 
