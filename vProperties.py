class Step1:

    '''Instantiation, from vProperties import Step1, x = Step1("virus"), x.s'''

    def __init__(self, species):
        if species == "virus":
            self.s = "virus"
        elif species == "bacteria":
            self.s= "bacteria"

class Step2:
    def __init__(self, species):
        if species == "virus":
            self.s = "virus"
        elif species == "bacteria":
            self.s= "bacteria"