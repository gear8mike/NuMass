def bubble_sort(array):
    n = len(array)
    cond = True
    while cond:
        cond = False
        for i in range(n-1):
            if array[i] > array[i+1]:
                array[i], array[i+1] = array[i+1], array[i]
                cond = True  # If a swap occurred, set condition to True to continue the loop

if __name__ == '__main__':
    array = [9, 5, 1, -8, 3, 7, 0, 25, 2]
    print(f'Unsorted array: {array}')
    bubble_sort(array)
    print(f'Sorted array: {array}')