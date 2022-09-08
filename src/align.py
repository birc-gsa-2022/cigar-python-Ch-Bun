"""A module for translating between alignments and edits sequences."""


def get_edits(p: str, q: str) -> tuple[str, str, str]:
    """Extract the edit operations from a pairwise alignment.

    Args:
        p (str): The first row in the pairwise alignment.
        q (str): The second row in the pairwise alignment.

    Returns:
        str: The list of edit operations as a string.

    >>> get_edits('ACCACAGT-CATA', 'A-CAGAGTACAAA')
    ('ACCACAGTCATA', 'ACAGAGTACAAA', 'MDMMMMMMIMMMM')

    """
    assert len(p) == len(q)
    # FIXME: do the actual calculations here
    #1. If the two strings are empty, you are done.
    if len(p) == 0:
        return('', '', '')
    # Otherwise:
    cigar = ''
    for i in range(len(p)):
        #2. If the two strings both have non-gaps at the front, add the letters there to the
        #corresponding output and put an `M` in the edits string.
        if (p[i] != '-' and q[i] != '-'):
            cigar = ''.join((cigar, 'M'))
        #3. If the first string has a gap, then the other doesn't (that is an invariant we will
        #insist on), so you add the second's string's letter to its corresponding output, and you
        #add an `I` to the edits.
        elif (p[i] == '-'):
            cigar = ''.join((cigar, 'I'))
        #4. If the second string has a gap, then add the first string's letter to its output and add
        #a 'D' to the edits.
        else:
            cigar = ''.join((cigar, 'D'))

    return (p.replace("-",""), q.replace("-",""), cigar)


def local_align(p: str, x: str, i: int, edits: str) -> tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> local_align("ACCACAGTCATA", "GTACAGAGTACAAA", 2, "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')
    """
    # FIXME: Compute the alignment rows
    x_upper_bound = edits.count('M') + edits.count('I')+i
    p_, x_ = align(p, x[i:x_upper_bound], edits)
    return p_ , x_


def align(p: str, q: str, edits: str) -> tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The first sequence to align.
        q (str): The second sequence to align
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> align("ACCACAGTCATAAA", "ACAGAGTACAAA", "MDMMMMMMIMMMMDD")
    ('ACCACAGT-CATAAA', 'A-CAGAGTACAAA--')

    """
    # FIXME: Compute the alignment rows

    for i in range(len(edits)):
    #3. If you have an `I` edit, you emit a `-` to the first string's output and the second string's letter
    #to the other output.
        if(edits[i] == 'I'):
            p = '-'.join((p[:i],p[i:]))   
    #4. If you have a `D` edit, you copy the first string's letter to its output and you emit a `-`
    #for the second string.
        elif (edits[i] == 'D'):
            q = '-'.join((q[:i],q[i:])) 
    return p, q


def edit_dist(p: str, x: str, i: int, edits: str) -> int:
    """Get the distance between p and the string that starts at x[i:]
    using the edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        int: The distance from p to x[i:?] described by edits

    >>> edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM")
    5
    """
    # FIXME: Compute the edit distance
    p_align, x_align = local_align(p, x, i, edits)
    edit_distance = 0
    for i in range(len(p_align)):
        if (p_align[i] != x_align[i]):
            edit_distance += 1
    return edit_distance

#Test
#def main():
#    print(get_edits('ACCACAGT-CATA', 'A-CAGAGTACAAA'))
#    print(get_edits(align("ACCACAGTCATA", "ACAGAGTACAAA", "MDMMMMMMIMMMM")))

#    print(align("ACCACAGTCATA", "ACAGAGTACAAA", "MDMMMMMMIMMMM"))
#    print(align("ACCACAGTCATAAA", "ACAGAGTACAAA", "MDMMMMMMIMMMMDD"))
#    print(align("", "", ""))

#    print(local_align("ACCACAGTCATA", "GTACAGAGTACAAA", 2, "MDMMMMMMIMMMM"))
#    print(local_align("accaaagta","gtacaaatgtcca",2,"MDMMIMMMMIIM"))
    
#    print(edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM"))

#if __name__ == '__main__':
#    main()