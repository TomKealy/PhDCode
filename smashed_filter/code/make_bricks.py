def make_bricks(small, big, goal):
    if goal > big*5 + small:
        return False
    needed = goal%5
    if needed > small:
        return False
    return True
