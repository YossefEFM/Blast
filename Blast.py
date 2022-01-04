import pandas as pd

protiens = ['A','R','N','D','C','Q','E',
            'G','H','I','L','K','M','F',
            'P','S','T','W','Y','V','B','Z','X','*']

blosum62 = pd.read_csv("BLOSUM62.csv")
blosum62 = blosum62.rename(index= lambda s: protiens[int(s)])
blosum62_map = {}

for i in protiens:
    for j in protiens:
        blosum62_map[(i, j)] = blosum62.at[i, j]
del blosum62

def get_words(sequence, wordLen):
    words=[]
    for i in range(0, len(sequence) - wordLen - 1):
        word = sequence[i: i+wordLen]
        words.append((word, i))
    return words

def get_neighbors(words, wordLen, threshold):
    seeds=[]
    for (word, idx) in words:
        temp_w, len_w = word, len(word)
        for i in range(len_w):
            for aa in protiens:
                 temp_w = temp_w[:i] + aa + temp_w[i+1:]
                 score = 0
                 for j in range(wordLen):
                     score = score + blosum62_map[(word[j], temp_w[j])]
                 if score >= threshold:
                     seeds.append((temp_w, word, score, idx))
    return seeds

def get_hits(protien, seeds, wordLen):
    hits=[]
    for i in range(len(seeds)):
        (temp_w, _, _, _) = seeds[i]
        for j in range(len(protien)):
            if protien[j:j+wordLen] == temp_w:
                hits.append((i, j))
    return hits

def hspExtend(hits, protien, wordLen):
    hsp = []
    for i in range(len(hits)):
        seedIdx = hits[i][0]

        leftExt = hits[i][1]
        rightExt = hits[i][1]+wordLen-1
        leftQuery = seeds[seedIdx][3]
        rightQuery = seeds[seedIdx][3]+wordLen-1

        score = seeds[seedIdx][2]
        scoreX = score

        while(scoreX - score < threshold):
            if(leftExt >= 0 and leftQuery >= 0):
                score += blosum62_map[(protien[leftExt], sequence[leftQuery])]
                leftQuery-=1
                leftExt-=1
            if(rightExt < len(protien) and rightQuery < len(sequence)):
                score += blosum62_map[(protien[rightExt], sequence[rightQuery])]
                rightExt+=1
                rightQuery+=1
            if leftExt < 0 or leftQuery < 0:
                if rightExt >= len(protien) or rightQuery >= len(sequence):
                    break
            if(scoreX < score):
                scoreX = score
        hsp.append([rightExt-1, leftExt+1, rightQuery-1, leftQuery+1, score])
    return hsp

def overLap(hsp, protien):
    NewHSP=[]
    idx=0
    while(idx<len(hsp)):
        frExt=hsp[idx][0]
        srExt=hsp[idx+1][0]

        flExt=hsp[idx][1]
        slExt=hsp[idx+1][1]

        if frExt>=slExt:
            frQuery=hsp[idx][2]
            srQuery=hsp[idx+1][2]

            flQuery=hsp[idx][3]
            slQuery=hsp[idx+1][3]

            if frQuery >= slQuery:
                if flQuery == slQuery and frQuery == srQuery:
                    NewHSP.append(hsp[idx])
                else:
                    if flQuery == slQuery:
                        if frQuery > srQuery:
                            NewHSP.append(hsp[idx])
                        else:
                            NewHSP.append(hsp[idx + 1])
                    else:
                        nSc = hsp[idx][4] + hsp[idx + 1][4]
                        counter = slExt
                        for n in range(slQuery, frQuery):
                            nSc -= blosum62_map[protien[counter]][sequence[n]]
                            counter += 1
                        NewHSP.append([srExt, flExt, srQuery, flQuery, nSc])
                idx += 1
            else:
                NewHSP.append(hsp[idx])
        else:
            NewHSP.append(hsp[idx])
        idx += 1
    return NewHSP

def display(NewHSP, protien):
    for z in range(len(NewHSP)):
        x = NewHSP[z][0]
        y = NewHSP[z][1]

        print("Protien : ", protien [y : x+1])
        print("Query : ", sequence[NewHSP[z][3]:NewHSP[z][2]+1])
        print("Score : ", NewHSP[z][4])



sequence = 'PQGEFG'
protien = 'PYGPWGETPCGWFRRQGEHADGKRGEFQWEAAAAAPWG'

wordLen = 3
threshold = 10

words = get_words(sequence, wordLen)
seeds = get_neighbors(words, wordLen, threshold)

hits = get_hits(protien, seeds, wordLen)
HSP = hspExtend(hits, protien, wordLen)
final = overLap(HSP, protien)
display(final, protien)