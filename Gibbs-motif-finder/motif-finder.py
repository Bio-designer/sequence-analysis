import random
import numpy
import math
import time


##this program now is designed only for finding motifs in DNA sueqnce, the alphabet can be extended easily in future to support more choices

##alphabet
abet= {
    "A":0,
    "T":1,
    "C":2,
    "G":3
}
abet_array=["A","T","C","G"]

##set the frequence of every char
##set it equally or you can compute according to (#x_i in entrie input sequence)/length of sequence
freq={
    "A":0.25,
    "T":0.25,
    "C":0.25,
    "G":0.25
}


##iterations
ITERATIONS=100000
##IC threshold for every position
IC_threshold=1.6
##fixed length for the motif
ml=6






##compute the PWM for a given sequence set
def getPWM(sss):
    ##PWM: position weight matrices
    pwm = numpy.zeros((len(abet), ml))
    for i in range(0,len(sss)):
        j=0
        for char in sss[i]:
            pwm[abet[char]][j] += 1
            j=j+1
    ##pay attention, here is len(sss), not len(sss)-1
    pwm=pwm/len(sss)

    return pwm





##compute the score for a given sequence according to a pwm
#the candidate sequence is the one to choose a new starting position
def getScore(candidate_sequence,pwm):
    ##since j is not included in the range(i,j), so here is len()-ml+1, rather than len()-ml
    ##please pay attention
    candidate_seq = [candidate_sequence[j:j + ml] for j in range(0, len(candidate_sequence) - ml + 1)]

    ##store every position's probablity
    prob=[]


    for j in range(0, len(candidate_sequence) - ml + 1):

        ##initiate the score
        score = 1
        ##counting position
        pos = 0
        for char in candidate_sequence[j:j+ml]:
            temp = pwm[abet[char]][pos] / freq[char]
            # if temp=0, just ignore it
            if temp != 0:
                ##does the base for the log function matters? 2 or 4
                #score += math.log(temp, len(abet))
                score*=temp
                #print(score)

            # score += math.log(pssm[abet[char]][pos] / freq[char],4)
            pos += 1
        prob.append(score)
    normed_prob=list(prob/sum(prob))

    return normed_prob



##compute the information content for a given motif sequence set
##actually, it is Kullback-Leiber distance

def getIC(motif_seq):
    ##first compute the pwm for the motif_seq
    pwm_motif=getPWM(motif_seq)
    ##calculate the IC for a given PWM

    ##set IC to zero
    IC=0

    for i in range(0, len(pwm_motif)):

        pos = 0
        for j in range(0,ml):
            temp=pwm_motif[i][j]/freq[abet_array[i]]
            if temp!=0:
                IC+=pwm_motif[i][j]*math.log(temp,2)
                ##please compare IC with score
            pos+=1
    return IC


#print(motif(seq))


##smapling new starting position according to its probality distribution
##ref to:
## http://stackoverflow.com/questions/479236/how-do-i-simulate-biased-die-in-python
def getSite(massDist):
    randRoll=random.random() #in [0,1)
    sum=0
    ##starting from 0, not 1
    result=0
    for mass in massDist:
        sum +=mass
        if randRoll<sum:
            return result
        result+=1



##output a sequence represents the most likely motif
def getMotif(PWM):

    mtx=numpy.transpose(PWM)

    motif=[]
    IC_unit=[]

    for pos in range(0,ml):
        temp=list(mtx[pos])
        max_val=max(temp)

        temp_IC_unit=0
        ##for every position in the motif sequence sets, calculate the corresponding IC
        for i in range(0,len(temp)):
            ttt=temp[i]/freq[abet_array[i]]
            ##if zero, log function doesn't work
            if ttt!=0:
                temp_IC_unit+=temp[i]*math.log(ttt,2)

        IC_unit.append(temp_IC_unit)


        ##check if the max_val appears more than once
        duplicate_check=[]
        for i in range(0,len(temp)):
            if (temp[i]==max_val):
                duplicate_check.append(abet_array[i])

        if (len(duplicate_check)>1):
            ##return all the equal case, rather than return one choice randomly
            motif.append(''.join(duplicate_check))
        else:
            motif.append(abet_array[temp.index(max_val)])

    return motif,IC_unit



















if __name__ == "__main__":

    ##read sequences for finding motif from input file
    seq = []
    with open("~/sequence.fa") as sequence_input:
        for line in sequence_input:
            seq.append(line.replace('\n', ''))

    sequence_input.close()

    ##random starting positions
    ##because there is another random which is used much more frequently, making the use here meaningless
    random.seed(time)
    rsp = [random.randint(0, len(x) - ml) for x in seq]


    maxIC=0
    ##recore the best starting positions
    best_start_position= [0 for x in seq]
    IC_min=IC_threshold * ml


    ##do interations
    for ITER in range(0,ITERATIONS):
        ##generate a random number i to ignore sequence i

        i = random.randrange(0, len(seq))

        ##get the sequence segment except for seq i
        sss = [x[j:j + ml] for q, (x, j) in enumerate(zip(seq, rsp)) if q != i]
        ##compute the PWM for the sequence segment set
        pwm = getPWM(sss)

        ##the score distribution of every possible starting position of the chosen sequence
        score_prob = getScore(seq[i], pwm)
        #print(score_prob)

        ##update the ith position according to the score distribution
        ##randomly choose one according to the probablity distribution
        rsp[i]=getSite(score_prob)

        ##this is important
        ##or you can just choose the best
        #best = score_prob.index(max(score_prob))
        #rsp[i] = best

        ##generate the current motif sequence set
        mss=[x[j:j + ml] for q, (x, j) in enumerate(zip(seq, rsp)) ]
        current_IC=getIC(mss)

        if current_IC>maxIC:
            maxIC=current_IC
            best_start_position=list(rsp)
        #since one drawback of Kullback-Leibler has no boundary, so the below choice is not easy to satisfy every special case
        if (float(current_IC) > IC_min):
            break #exist the loop, already find the satisfactory answer

    IC=maxIC
    mss = [x[j:j + ml] for q, (x, j) in enumerate(zip(seq, best_start_position))]
    ##recored the final chosen motif sequence set and its PWM
    PWM=getPWM(mss)


    print("Final PWM:")
    print(PWM)
    print("Fianl motif sequence set:")
    print(mss)

    print("Total information content for the motif:")
    print(IC)

    print("Iteration times:")
    print(ITERATIONS)

    print("The starting position of the motifs in every sequence:")
    print(best_start_position)

    finalmotif, ICU=getMotif(PWM)
    print("The best matched motif:")
    print(finalmotif)
    print("Information content per site:")
    print(ICU)














