class VariableGenerator(object): 
     def __init__(self, prefix): 
         self.__prefix = prefix 
     @cached_method 
     def __getitem__(self, key): 
         return SR.var("%s%s"%(self.__prefix,key)) 

