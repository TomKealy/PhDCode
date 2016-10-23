class Tree:
    def __init__(self, x, left = None, right = None):
        self.val = x
        self.left = left
        self.right = right
        
    def Traverse(self):
        print self.val
        if self.left:
            self.left.Traverse()
        if self.right:
            
            self.right.Traverse()
            
    def depth(self):
        left_depth = self.left.depth() if self.left else 0
        right_depth = self.right.depth() if self.right else 0
        return max(left_depth, right_depth) + 1

tree = Tree('+', Tree(1), Tree('*', Tree(2), Tree(3)))
