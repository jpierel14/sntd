

class fits:
    '''
    This is going to be a fits object that will contain a dictionary
    of models with keys as model names. Each value will then be a fitted model
    class from sncosmo.

    If microlensing is added, I guess run through pycs first to get initial time
    delays and microlensing effects, write out the result, and run through sncosmo
    to get best fit model?
    '''