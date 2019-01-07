class Interval:
    '''
    Window in which we compute muational signature
    '''
    def __init__(self, start, end):
        
        self.start = start
        self.end = end
        
        
    def is_included(self, pos):
        
        if self.start <= pos <= self.end:
            return True
        return False