import time


def benchmark(func):
    def wrapper(*args, **kwargs):
        begin = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print('Execution time (' + func.__name__ + '): ' + str(end - begin) + 's')
        return result
    return wrapper

