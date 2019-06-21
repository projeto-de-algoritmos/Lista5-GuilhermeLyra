def minimum_penalty(x, y, p_mm, p_gap):
    m,n=len(x),len(y)
    dp=[[0 for i in range(n+m+1)] for j in range(n+m+1)]
    for i in range(n+m+1):
        dp[i][0]=i*p_gap;
        dp[0][i]=i*p_gap;
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            if(x[i-1]==y[j-1]):
                dp[i][j]=dp[i-1][j-1]
            else:
                dp[i][j]=min({dp[i - 1][j - 1] + p_mm, dp[i - 1][j] + p_gap, dp[i][j - 1] + p_gap})

    l = n+m
    i,j=m,n
    xpos=ypos=l
    xans=[0 for i in range(l+1)]
    yans=[0 for i in range(l+1)]
    while not(i==0 or j==0):
        if(x[i-1]==y[j-1]):
            xans[xpos]=ord(x[i-1])
            yans[ypos]=ord(y[j-1])
            i-=1
            j-=1
            xpos-=1
            ypos-=1
        elif (dp[i - 1][j - 1] + p_mm == dp[i][j]):
            xans[xpos] = ord(x[i - 1])
            yans[ypos] = ord(y[j - 1])
            i-=1
            j-=1
            xpos-=1
            ypos-=1
        elif (dp[i - 1][j] + p_gap == dp[i][j]):
            xans[xpos] = ord(x[i - 1])
            yans[ypos] = ord('_')
            i-=1
            xpos-=1
            ypos-=1
        elif (dp[i][j - 1] + p_gap == dp[i][j]):
            xans[xpos] = ord('_')
            yans[ypos] = ord(y[j - 1])
            j-=1
            xpos-=1
            ypos-=1
    while xpos > 0:
        if i > 0:
            xans[xpos]=ord(x[i])
            i-=1
        else:
            xans[xpos]=ord('_')
        xpos-=1

    while ypos > 0:
        if j > 0:
            yans[ypos]=ord(x[j])
            j-=1
        else:
            yans[ypos]=ord('_')
        ypos-=1

    idx = 1
    for i in range(l, 0, -1):
        if(chr(yans[i])=='_' and chr(xans[i])=='_'):
            idx+=1        
    
    print("menor penalidade: {0}".format(dp[m][n]))
    print(''.join([chr(a) for a in xans[idx:l+1]]))
    print(''.join([chr(a) for a in yans[idx:l+1]]))
    # for i in range(idx, l+1):
    #     print(chr(xans[i]), end='')
    # for i in range(idx, l+1):
    #     print(chr(yans[i]))
    
    return

# gene1 = "GCAAGGAGTCGACCGTCGCAAACAAACGCCATGTGTCTCATGTTCTCGTTCGGTCTATCCCAAAGCGGATCAGCTGGAAGCGGGATTCGAAAACAGCGACATCCCTGTCTTGGGCAAACGCATTGCACACATTCCTTCGGGTGAAACAATGAAGTTACACACTCAGTAAGGCTTCTTGGCTTATCACTTGTCATTGGGTGCTCAGTGGCCGATAGATGCAGTCTGTATATAAAGATGATCAACACAACCTATAATGACGGCCATGTAACACCCCTAAGCAAACACCCTGAACCAAAACATCCTCTACACTCTATCTAGTCAATCTTCGGAATGGTGAACACCTGCACCTATCTCCCCC";
# gene2 = "ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGTGTGCTGAGCCGAATGAAAGTTGCCTAATTACTAAGGTGTAGTTCCAGCATACCATACACCCTAACTCATACTACGGTAGGTAGATCTACTTACCTATGAACCTATATTGGTAGGTAGGTGAATATAAAATACAGCATGGAACATGTTTTTCATTAGCTGGTCTCTCATTCGTCCTTGTCCTAGGCCTTAAGGAATCCAGTATATGAAATAATCCCTCTTATCCATTTTCCTCCTATTCTTTTTCATTTCCCTCATCACTGCAACTCTAATCCTCGGGCTCACCCTCCCTGTGTCTCCTCGAAATGGTG";

# minimum_penalty(gene1, gene2, 1, 2)