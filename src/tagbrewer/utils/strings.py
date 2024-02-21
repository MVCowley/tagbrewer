def sliceIterator(lst, sliceLen):
    """
    Source: https://stackoverflow.com/questions/1335392/iteration-over-list-slices
    """
    for i in range(len(lst) - sliceLen + 1):
        yield lst[i:i + sliceLen]