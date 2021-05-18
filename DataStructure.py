from dataclasses import dataclass

@dataclass
class Species(): 
    mass:float = None
    name:str = None
    area:float = None
    intensity:float = None
    timepoint:list = None

        
class TimePoint():
    def __init__(self):
        self.species = {}

class Reaction():
    def __init__(self, ID, species = []):
        self.ID = ID
        self.time_series = {}
        self.metadata = {}
        
# class MetaData():
#     def __init__(self):
#         pass
    
