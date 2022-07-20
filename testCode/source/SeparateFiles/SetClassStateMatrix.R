#### State Objects (matrix that holds the population state) ####

#Matrix that holds the population state in terms of infected, susceptible, recovered, etc.
setClass("ID.state.matrix",
         slots = list(n.epi.class = "numeric", #number of different
                      #epi classes
                      epi.class = "numeric", #the epid class of each row,
                      #between 1 and n.epi.class
                      epi.class.label = "character", #label for each epi class
                      n.age.class = "numeric"#the number of age classes
         ),
         contains="matrix")
