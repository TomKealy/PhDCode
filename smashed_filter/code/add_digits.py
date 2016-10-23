def addDigits(num):
        """
        :type num: int
        :rtype: int
        """
        s = len(str(num))
        next = sum(map(int, str(num)))
        while s != 1:
            next = sum(map(int, str(next)))
            s = len(str(next))
        return next
