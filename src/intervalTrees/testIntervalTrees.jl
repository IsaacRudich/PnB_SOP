#include("treesMain.jl")
include("../utilities.jl")


function testIntervalTrees()
    #test add nodes
    rootNode = IntervalTreeNode(5,10,1)
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(15,25))
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(1,12,2))
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(8,16))
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(14,20))
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(18,21))
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(2,8))
    addNodeToIntervalTree!(rootNode, IntervalTreeNode(2.1,8))

    print("Ordered Intervals From Tree \n")
    printTree(rootNode)
    print("\n\n")

    #testTreeCopy, display at end
    treeCopy = deepCopyIntervalTreeNodes(rootNode, false)
    treeCopyMerged = deepCopyIntervalTreeNodes(rootNode, true)


    #test query interval
    print("Query Interval (8,10,1) \n")
    query1 = Interval(8,10)
    results = IntervalTreeNode[]
    findIntersectingIntervals!(query1, rootNode, results)
    for result in results
        print(result,"\n")
    end

    #test query interval again
    query2 = Interval(20,22)
    results = IntervalTreeNode[]
    print("\n\n")
    print("Query Interval (20,22,1) \n")
    findIntersectingIntervals!(query2, rootNode, results)
    for result in results
        print(result,"\n")
    end
    print("\n\n")



    #test delete interval
    #=
    query3 = Interval(5,5,1)
    results = IntervalTreeNode[]
    findIntersectingIntervals!(query3, rootNode, results)
    print("Tree Without Intervals Intersecting '[]' \n")
    for result in results
        print("Deleting: ", result, "\n")
        rootNode = deleteNodeFromIntervalTree!(result)
    end
    printTree(rootNode)
    =#


    #test split on point
    print("Tree Without 8 in it \n")
    rootNode = removePointFromIntervalTree!(rootNode, 8)
    print("Ordered Intervals From Tree \n")
    printTree(rootNode)
    print("\n\n")

    #test split on point
    print("Tree Without [5,10] in it \n")
    rootNode = removeIntervalFromIntervalTree!(rootNode, Interval(5,10))
    print("Ordered Intervals From Tree \n")
    printTree(rootNode)
    print("\n\n")

    print("Impose max of 15 \n")
    rootNode = imposeMaximumOnIntervalTree!(rootNode, 15)
    print("Ordered Intervals From Tree \n")
    printTree(rootNode)
    print("\n\n")

    print("Impose min of 2 \n")
    #rootNode = imposeMinimumOnIntervalTree!(rootNode, 2)
    print("Ordered Intervals From Tree \n")
    printTree(rootNode)
    print("\n\n")

    print("Impose open max of 14 \n")
    rootNode = imposeOpenMaximumOnIntervalTree!(rootNode, 14)
    print("Ordered Intervals From Tree \n")
    printTree(rootNode)
    print("\n\n")

    println("Impose open min of 3")
    rootNode = imposeOpenMinimumOnIntervalTree!(rootNode, 3)
    println("Ordered Intervals From Tree")
    printTree(rootNode)
    print("\n\n")

    println("Add the interval [8,15] (1) with merging")
    rootNode = addIntervalToIntervalTreeWithMerging!(rootNode, Interval(8,15))
    println("Ordered Intervals From Tree")
    printTree(rootNode)
    print("\n\n")

    println("Add the interval [26,30] (1) with merging")
    rootNode = addIntervalToIntervalTreeWithMerging!(rootNode, Interval(26,30))
    println("This is the 'Post Manipulation' Tree")
    printTree(rootNode)
    print("\n\n")

    println("Original Tree Copy")
    println("Ordered Intervals From Tree")
    printTree(treeCopy)
    print("\n\n")

    println("Original Tree Merged")
    println("Ordered Intervals From Tree")
    printTree(treeCopyMerged)
    print("\n\n")

    println("Union of Original and Post-Manipulation")
    println("Ordered Intervals From Tree")
    union = unionIntervalTrees(treeCopy, rootNode, false)
    printTree(union)
    print("\n\n")

    println("Union of Original and Post-Manipulation With Merging")
    println("Ordered Intervals From Tree")
    union2 = unionIntervalTrees(treeCopy, rootNode, true)
    printTree(union2)
    print("\n\n")

    ex1 = IntervalTreeNode(26,30,1)
    addNodeToIntervalTree!(ex1, IntervalTreeNode(15,25))
    addNodeToIntervalTree!(ex1, IntervalTreeNode(100,200,3))
    addNodeToIntervalTree!(ex1, IntervalTreeNode(300,400,6))

    ex2 = IntervalTreeNode(5,10,1)
    addNodeToIntervalTree!(ex2, IntervalTreeNode(15,25))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(1,12,2))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(8,16))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(14,20))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(18,21))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(2,8))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(90.4,150,1.6))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(169,200,0.6))
    addNodeToIntervalTree!(ex2, IntervalTreeNode(310,400,10))

    println("Ordered Intervals From New Tree 1")
    printTree(ex1)
    print("\n\n")
    println("Ordered Intervals From New Tree 2")
    printTree(ex2)
    print("\n\n")


    println("Intersection of Previous 2 Trees")
    println("Ordered Intervals From Tree")
    intersection = intersectIntervalTrees(ex1, ex2)
    printTree(intersection)
    print("\n\n")

    println("Test of creating a tree from an array")
    a = [1,2,3,4,5,7,8,9,10]
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of creating a tree from an array")
    a = [1,2,3,4,5,7,8,9,11]
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of creating a tree from an array")
    a = [1,2,5,7,4,11,11.5,13.4,13.2]
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of creating a tree from an array")
    a = [1,1,1,1]
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of creating a tree from an array")
    a = [1,10]
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of creating a tree from an array")
    a = [1]
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of creating a tree from a set")
    a = Set((1,2,3,4,5,7,8,9,10))
    println(a)
    treeFromSet = setToIntervalTree(a)
    printTree(treeFromSet)
    print("\n\n")

    println("Test of disjoint union")
    b = Set((2,4,8,10,12))
    println(a)
    println(b)
    println("Expected Answer: {12} {1,3,5,7,9}")
    treeFromSet2 = setToIntervalTree(b)
    setA, setB = disjunctiveUnion(treeFromSet, treeFromSet2)
    println(setA)
    println(setB)
end
