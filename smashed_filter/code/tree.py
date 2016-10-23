class TreeNode:
    """Simple Tree class with arbirary number of children"""
    def __init__(self, data):
        self.data = data
        self.children = []

    def add_children(self, child):
        self.children.append(child)

    def __str__(self):
        return str(self.data)

    def print_tree(self, root):
        if root is None:
            return
        print root.data
        for child in root.children:
            self.print_tree(child)

r = TreeNode(0)
ch1 = TreeNode(1)
ch2 = TreeNode(2)
ch3 = TreeNode(3)
r.add_children(ch1)
r.add_children(ch2)
r.add_children(ch3)

ch4 = TreeNode(4)
ch1.add_children(ch4)
