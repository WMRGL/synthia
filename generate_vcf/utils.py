import multiprocessing


class MultiprocessCounter(object):
    def __init__(self, initval=0):
        self.counter = multiprocessing.RawValue('i', initval)
        self.lock = multiprocessing.Lock()

    def increment(self, val=1):
        with self.lock:
            self.counter.value += val

    @property
    def value(self):
        return self.counter.value
