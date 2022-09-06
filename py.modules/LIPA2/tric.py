#

class Tricky:
    def __init__(self,x):
        self.data=x

    def __setstate__(self, d):
        print('setstate happening')
        self.data = 10

    def __getstate__(self):
        return self.data
        print('getstate happening')
