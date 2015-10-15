import profile
import itertools
import collections
import string
import re
import sys
from collections import *

# megkapjuk a beolvasott file osszes token-et (szavat) kozpontozas nelkul
def get_tokens_list(file_name):
    tokens_list = []
    punct = set(string.punctuation)
    '!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890'

    letters = set('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890]')

    with open(file_name,'r') as f:
         for line in f:
            for word in line.split():
                word = re.sub('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890]','',word)
                if not word.startswith('-') and word != '' and not (letters & set(word)):
                    word = word.lower()
                    #word = unicode(word,'utf-8')
                    tokens_list.append(word+' ')
    return tokens_list


# megkapjuk az osszes betu szamat
def get_all_letters_nr(tokens_list):
    all_letters = 0

    for word in tokens_list:
        lenght = len(word)
        all_letters += lenght

    return all_letters

# osszeszamolja, hogy j egymast koveto karakterbol mennyi van
# a szovegben
def n_gram_letter_counter(j,tokens_list):
    ngram_list = collections.defaultdict(int)

    for word in tokens_list:
        lenght = len(word)
        i = 0
        while lenght-j+1 > i:
            ngram_list[word[i:i+j]] += 1
            i += 1

    return ngram_list

################ UNSMOOTHING probabilities ######################

def ngram_unsmoothing_prob(i,ngram_letter_counter,n_1gram_letter_counter,all_letters_nr):
    ngram_prob_list = collections.defaultdict(int)

    if i == 1:
        for ngram in ngram_letter_counter:
            ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0) / all_letters_nr
        return ngram_prob_list
    else:
        for ngram in ngram_letter_counter:
            if n_1gram_letter_counter[ngram[0:i-1]] != 0:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0) / n_1gram_letter_counter[ngram[0:i-1]]
            else:
                ngram_prob_list[ngram] = 0.0
        return ngram_prob_list
######################################################################

#################### ADD-ONE probabilities ############################

def ngram_add_one_prob(i,ngram_letter_counter,n_1gram_letter_counter,all_letters_nr,letters_types):
    ngram_prob_list = collections.defaultdict(int)

    if i == 1:
        for ngram in ngram_letter_counter:
            ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+1.0) / (all_letters_nr+letters_types)
        return ngram_prob_list
    else:
        for ngram in ngram_letter_counter:
            ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+1.0) / (n_1gram_letter_counter[ngram[0:i-1]]+letters_types)
        return ngram_prob_list
########################################################################

############## WITTEN-BELL DISCOUNTING probabilities ####################

def ngram_witten_bell_prob(i,ngram_letter_counter,n_1gram_letter_counter):
    n_1gram_z_list = collections.defaultdict(int)
    n_1gram_watched_types = collections.defaultdict(int)
    tokens_nr = 0
    ngram_prob_list = collections.defaultdict(int)
    z = 0

    if i != 1:
        for ngram in ngram_letter_counter:
            if ngram_letter_counter[ngram] == 0:
                n_1gram_z_list[ngram[0:i-1]] += 1

        for n_1gram in n_1gram_letter_counter:
            n_1gram_watched_types[n_1gram] = len(n_1gram_letter_counter) - n_1gram_z_list[n_1gram]
            tokens_nr = tokens_nr + n_1gram_letter_counter[n_1gram]
    else:
        for ngram in ngram_letter_counter:
            if(ngram_letter_counter[ngram] == 0):
                z += 1
            else:
                tokens_nr = tokens_nr + ngram_letter_counter[ngram]

        watched_types_nr = len(ngram_letter_counter) - z

    for ngram in ngram_letter_counter:
        if i != 1:
            if(ngram_letter_counter[ngram] == 0):
                ngram_prob_list[ngram] = (n_1gram_watched_types[ngram[0:i-1]]+0.0)/(n_1gram_z_list[ngram[0:i-1]]*(tokens_nr+n_1gram_watched_types[ngram[0:i-1]]))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0)/(n_1gram_letter_counter[ngram[0:i-1]]+n_1gram_watched_types[ngram[0:i-1]])
        else:
            if(ngram_letter_counter[ngram] == 0):
                ngram_prob_list[ngram] = (watched_types_nr + 0.0) / (z * (tokens_nr + watched_types_nr))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram] + 0.0) / (tokens_nr + watched_types_nr)

    return ngram_prob_list
#############################################################################################

############################ GOOD-TURING DISCOUNTING ##################################

# az occuring listbe szamolja, hogy az hany darab ngram van x elofordulassal
def ngram_occuring_list(ngram_letter_counter,k):
    occuring_list = collections.defaultdict(int)

    for x in range(0,k):
        for ngram in ngram_letter_counter:
            if(ngram_letter_counter[ngram] == x):
                occuring_list[x] += 1
    return occuring_list

def ngram_good_turing_discounting(ngram_occuring_list,k):
    smoothed_count_c_list = collections.defaultdict(int)

    for x in range(0,k-1):
        if ngram_occuring_list[x] != 0:
            smoothed_count_c_list[x] = (x+1.00)*ngram_occuring_list[x+1]/ngram_occuring_list[x]
        else:
            smoothed_count_c_list[x] = 0.0
    return smoothed_count_c_list
########################################################################################

######################################## BACKOFF##########################################

def prob_tilde(ngram_letter_counter,n_1gram_letter_counter,all_letters_nr):
    prob_tilde_list = collections.defaultdict(int)
    count_list = collections.defaultdict(int)
    tokens_nr = 0
    z = 0

    for ngram in ngram_letter_counter:
        if(ngram_letter_counter[ngram] == 0):
            z += 1
        else:
            tokens_nr = tokens_nr + ngram_letter_counter[ngram]

        watched_types_nr = len(ngram_letter_counter) - z

    for ngram in ngram_letter_counter:
        if ngram_letter_counter[ngram] == 0:
            count_list[ngram] = watched_types_nr*tokens_nr/(z*(tokens_nr+watched_types_nr+0.00))
        else:
            count_list[ngram] = ngram_letter_counter[ngram]*tokens_nr/(tokens_nr+watched_types_nr+0.00)

    for ngram in ngram_letter_counter:
        if len(ngram) != 1:
            if n_1gram_letter_counter[ngram[0:len(ngram)-1]] != 0:
                prob_tilde_list[ngram] = (count_list[ngram]+0.0)/n_1gram_letter_counter[ngram[0:len(ngram)-1]]
            else:
                prob_tilde_list[ngram] = 0.0
        else:
            prob_tilde_list[ngram] = (count_list[ngram]+0.0)/all_letters_nr

    return prob_tilde_list

def alpha(prob_tilde_list,prob_1tilde_list,ngram_letter_counter):
    sum_n_1gram_list = 0
    sum_ngram_list = 0

    for ngram in ngram_letter_counter:
        if prob_tilde_list[ngram] > 0:
            sum_ngram_list = sum_ngram_list + prob_tilde_list[ngram]
            sum_n_1gram_list = sum_n_1gram_list + prob_1tilde_list[ngram[1:]]

    alpha = (1.0-sum_ngram_list)/(1-sum_n_1gram_list)

    return alpha

def ngram_backoff(ngram_letter_counter,prob_tilde,n_1prob_tilde,n_2prob_tilde,n_1gram_letter_counter,alpha,alpha_1):
    backoff_prob_list = collections.defaultdict(int)

    for ngram in ngram_letter_counter:
        if ngram_letter_counter[ngram] > 0:
            backoff_prob_list[ngram] = prob_tilde[ngram]
        else:
            if n_1gram_letter_counter[ngram[1:]] > 0:
                backoff_prob_list[ngram] = alpha * n_1prob_tilde[ngram[1:]]
            else:
                backoff_prob_list[ngram] = alpha_1 * n_2prob_tilde[ngram[2:]]

    return backoff_prob_list

##############################################################################################

####################if the letter occures less than 30 the letter will be *##################
def changeSmallNrOfCharacters(ngram_letter_counter):
    i =0
    new_ngram_letter_counter=collections.defaultdict(int)
    for letter in ngram_letter_counter:
        if ngram_letter_counter[letter] <= 30:
            i=i+ngram_letter_counter[letter]
            ngram_letter_counter[letter]=0
    ngram_letter_counter['*']=i

    for letter in ngram_letter_counter:
        if ngram_letter_counter[letter] > 0:
            new_ngram_letter_counter[letter] = ngram_letter_counter[letter]

    return new_ngram_letter_counter

########################### TRAIN #############################
def replace_file_without_special_char(filename, unigram_letter_nr):
    with open(filename,"a") as myfile:
        for word in tokens_list:
            myfile.write(word)

def train_char_ngram(fname,order,orderplus1_ngram_prob):
    data = file(fname).read()
    lm = defaultdict(Counter)
    i=0

    for i in xrange(len(data)-order):
        history,char =data[i:i+order],data[i+order]
        if char != " " and not (history.startswith(" ") or history.endswith(" ")):
                word=history+char
                lm[history][char]=orderplus1_ngram_prob[word]
    return lm

def print_train(train_list):
    sorted(train_list)
    for letter in train_list:
        print letter
        print str(train_list[letter].most_common(10))

##################### FOR PRINTING #####################################

def print_unigram_and_probs(sorted_letter_counter,ngram_unsmoothing_prob,ngram_add_one_prob,ngram_witten_bell):
    for value,key in sorted_letter_counter:
        value.encode('utf-8')
        print 'unigram: {0}-{1}'.format(value,key)
        print 'unsmoothing prob: ' + str(ngram_unsmoothing_prob[value])
        print 'add-one smoothing prob:' + str(ngram_add_one_prob[value])
        print 'witten-bell discounting prob: ' + str(ngram_witten_bell[value])

############################### MAIN ############################
if __name__ == '__main__':
    file_name = sys.argv[1]

    tokens_list = get_tokens_list(file_name)
    #print sorted(tokens_list)
    all_letters_nr = get_all_letters_nr(tokens_list)
    #print all_letters_nr
    k = 11

    ############## UNIGRAM #################
    unigram_letter_counter = n_gram_letter_counter(1,tokens_list)
    letters_type = len(unigram_letter_counter)

    #sorted_unigram_letter_counter = sorted(unigram_letter_counter.iteritems(), key=lambda (k,v):v,reverse=True)

    #print sorted_unigram_letter_counter
    #print 'unigram letter counter: ' + str(sorted(unigram_letter_counter.iteritems(),key=lambda (k,v): v,reverse=True))

    #print 'unigram letter counter: ' + str(sorted(changeSmallNrOfCharacters(unigram_letter_counter).iteritems(),key=lambda (k,v): v,reverse=True))

    #unigram_unsmoothing_prob = ngram_unsmoothing_prob(1,unigram_letter_counter,unigram_letter_counter,all_letters_nr)
    #print 'unigram unsmoothing probabilities: ' + str(sorted(unigram_unsmoothing_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    #unigram_add_one_prob = ngram_add_one_prob(1,unigram_letter_counter,unigram_letter_counter,all_letters_nr,letters_type)
    #print 'unigram add-one probabilities: ' + str(sorted(unigram_add_one_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    #unigram_witten_bell = ngram_witten_bell_prob(1,unigram_letter_counter,1)
    #print 'unigram witten-bell discounting: ' + str(sorted(unigram_witten_bell.iteritems(),key=lambda (k,v): v,reverse=True))

    ############## BIGRAM ##################
    #bigram_letter_counter = n_gram_letter_counter(2,tokens_list)
    #print bigram_letter_counter
    #print 'bigram letter counter: ' + str(sorted(bigram_letter_counter.iteritems(),key=lambda (k,v): v,reverse=True))

    #bigram_unsmoothing_prob = ngram_unsmoothing_prob(2,bigram_letter_counter,unigram_letter_counter,all_letters_nr)
    #print 'bigram unsmoothing probabilities: ' + str(sorted(bigram_unsmoothing_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    #bigram_add_one_prob = ngram_add_one_prob(2,bigram_letter_counter,unigram_letter_counter,all_letters_nr,letters_type)
    #print 'bigram add-one probabilities: ' + str(sorted(bigram_add_one_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    #bigram_witten_bell = ngram_witten_bell_prob(2,bigram_letter_counter,unigram_letter_counter)
    #print 'bigram witten-bell discounting: ' + str(sorted(bigram_witten_bell.iteritems(),key=lambda (k,v): v,reverse=True))

    ############## TRIGRAM #################
    trigram_letter_counter = n_gram_letter_counter(3,tokens_list)
    #print trigram_letter_counter
    #print 'trigram letter counter: ' + str(sorted(trigram_letter_counter.iteritems(),key=lambda (k,v): v,reverse=True))

    #trigram_unsmoothing_prob = ngram_unsmoothing_prob(3,trigram_letter_counter,bigram_letter_counter,all_letters_nr)
    #print 'trigram unsmoothing probabilities: ' + str(trigram_unsmoothing_prob)

    #trigram_add_one_prob = ngram_add_one_prob(3,trigram_letter_counter,unigram_letter_counter,all_letters_nr,letters_type)
    #print 'trigram add-one probabilities: ' + str(trigram_add_one_prob)

    #trigram_witten_bell = ngram_witten_bell_prob(3,trigram_letter_counter,bigram_letter_counter)
    #print 'trigram witten-bell discounting' + str(trigram_witten_bell)
    
    ########### TRAIN ###################################
    replace_file_without_special_char("train_without.txt",tokens_list)
    train_bigram_unsmoothing = train_char_ngram("train_without.txt",2,trigram_unsmoothing_prob)
    print_train(train_bigram_unsmoothing)