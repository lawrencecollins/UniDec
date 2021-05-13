from dataclasses import dataclass

@dataclass
class Species(): 
    mass:float = None
    name:str = None
    area:float = None
    intensity:float = None
        
class TimePoint():
    def __init__(self):
        self.species = {}

class Reaction():
    def __init__(self, ID):
        self.ID = ID
        self.time_series = {}
        self.metadata = {}
        
# class MetaData():
#     def __init__(self):
#         pass
    
