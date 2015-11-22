import profile
import itertools
import collections
import string
import re
import sys
from collections import *
from scipy.stats import rankdata
import os
import math
import pycld2 as cld2
from itertools import tee,islice

# megkapjuk a beolvasott file osszes token-et (szavat) kozpontozas nelkul
def get_tokens_list(i,input_file_name):
    tokens_list = collections.defaultdict(int)
    punct = set(string.punctuation)
    '!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890'

    letters = set('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890]')

    with open(input_file_name,'r') as infile:
         for line in infile:
            for word in line.split():
                if not word.startswith('-') and word != '':
                    word = re.sub('[!"#$%&\'()+,-./:;<=>?@[\\]^_`{|}~1234567890]','',word)
                    word = word.lower()
                    #word = unicode(word,'utf-8')
                    space=' '*i
                    if word != '':
                        tokens_list[space+word+space] += 1
    return tokens_list

# osszeszamolja, hogy j egymast koveto karakterbol mennyi van
# a szovegben
def uni_gram_letter_counter(j,tokens_list):
    ngram_list = collections.defaultdict(int)

    for word in tokens_list:
        lenght = len(word)
        i = 0
        while lenght-j+1 > i:
            ngram_list[word[i:i+j]] += tokens_list[word]
            i += 1
    alma = ngram_list[' ']-2/2
    ngram_list[' ']=alma

    return ngram_list

def ngram_letter_counter(j,unigram_letter_counter,tokens_list):
    vocabulary=[]
    ngram_letter_counter = collections.defaultdict(int)
    
    for unigram,value in unigram_letter_counter.iteritems():
        vocabulary.append(unigram);
    
    ngram_counter_list=collections.defaultdict(int)
    ngram_l=(tuple("".join(item) for item in itertools.product(vocabulary,repeat=j)))
    
    for ngram in ngram_l:
        ngram_letter_counter[ngram]=0
    
    for word,value in tokens_list.iteritems():
        lenght = len(word)
        i = 0
        while lenght-j+1 > i:
            ngram_letter_counter[word[i:i+j]] += value
            i += 1
            
    return ngram_letter_counter

def ngram_types(letter_counter):
    ngram_letter_types = collections.defaultdict(int)
    
    for ngram,value in letter_counter.iteritems():
        if value>0:
            ngram_letter_types[ngram]=value
    
    return ngram_letter_types
    

################ UNSMOOTHING probabilities ######################

def ngram_unsmoothing_prob(i,ngram_letter_counter,n_1gram_letter_counter,tokens_nr):
    ngram_prob_list = collections.defaultdict(int)

    if i == 1:
        for ngram,value in ngram_letter_counter.iteritems():
            ngram_prob_list[ngram] = (value+0.0) / tokens_nr
        return ngram_prob_list
    else:
        for ngram,value in ngram_letter_counter.iteritems():
            #print ngram
            #print ngram[0:i-1]
            if n_1gram_letter_counter[ngram[0:i-1]] != 0:
                #print value
                #print n_1gram_letter_counter[ngram[0:i-1]]
                ngram_prob_list[ngram] = (value+0.0) / n_1gram_letter_counter[ngram[0:i-1]]
                #print ngram_prob_list[ngram]
            else:
                ngram_prob_list[ngram] = 0.0
            #print '~'*80
        
    return ngram_prob_list
######################################################################

#################### ADD-ONE probabilities ############################

def ngram_add_one_prob(i,ngram_letter_counter,n_1gram_letter_counter,tokens_nr,letter_types_nr):
    ngram_prob_list = collections.defaultdict(int)
    
    for ngram,value in ngram_letter_counter.iteritems():
        if i == 1:
            ngram_prob_list[ngram] = (value+1.0) / (tokens_nr+letter_types_nr)
        else:
            #print ngram
            #print ngram[0:i-1]
            #print ngram_letter_counter[ngram]
            #print n_1gram_letter_counter[ngram[0:i-1]]
            #print letter_types_nrs
            ngram_prob_list[ngram] = (value+1.0) / (n_1gram_letter_counter[ngram[0:i-1]]+letter_types_nr)
        
    return ngram_prob_list
########################################################################

############## WITTEN-BELL DISCOUNTING probabilities ####################

def ngram_witten_bell_prob(i,ngram_letter_counter,tokens_nr,letter_types):
    z_list = collections.defaultdict(int)
    observed_type = collections.defaultdict(int)
    ngram_prob_list = collections.defaultdict(int)
    new_bigram_list = collections.defaultdict(int)
    z = 0
    
    for ngram,value in ngram_letter_counter.iteritems():
        if value != 0:
            new_bigram_list[ngram[0:i-1]] += 1
            
    for ngram,value in ngram_letter_counter.iteritems():
        if i!= 1:
            if value == 0:
                observed_type[ngram] = new_bigram_list[ngram[0:i-1]]
                z_list[ngram]=letter_types-new_bigram_list[ngram[0:i-1]]
        else:
            if(value == 0):
                z += 1

    for ngram,value in ngram_letter_counter.iteritems():
        if i != 1:
            if(value == 0):
                ngram_prob_list[ngram] = (observed_type[ngram]+0.0)/(z_list[ngram]*(tokens_nr+observed_type[ngram]))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0)/(tokens_nr+observed_type[ngram])
        else:
            if(ngram_letter_counter[ngram] == 0):
                ngram_prob_list[ngram] = (letter_types + 0.0) / (z * (tokens_nr + letter_typesletter_types))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram] + 0.0) / (tokens_nr + letter_types)
        
    return ngram_prob_list

def c_star(ngram_letter_counter,ngram_witten_bell,tokens_nr):
    c_star = collections.defaultdict(int)
    
    for ngram,value in ngram_letter_counter.iteritems():
        c_star[ngram]=tokens_nr*ngram_witten_bell[ngram]
    return c_star
            
#############################################################################################

######################################## BACKOFF##########################################
def unigram_prob_tilde(ungram_letter_counter,tokens_nr):
    unigram_prob_tilde_list = collections.defaultdict(int)
    
    for ngram,value in unigram_letter_counter.items():
        unigram_prob_tilde_list[ngram]=value/(tokens_nr+0.0)
    
    return unigram_prob_tilde_list

def prob_tilde(ngram_letter_counter,n_1gram_letter_counter,c_star):
    prob_tilde_list = collections.defaultdict(int)
    count_list = collections.defaultdict(int)
    denominator_part=0

    for ngram, value in ngram_letter_counter.items():
        #print 'ngram and value= ' + str(ngram) +':'+str(value)
        #counter_part = ngram_witten_bell_prob[ngram]*tokens_nr
        counter_part = c_star[ngram]
        #print counter_part
        #print ngram
        #print ngram[0:len(ngram)-1]
        #print n_1gram_letter_counter[ngram[0:len(ngram)-1]]
        denominator_part = n_1gram_letter_counter[ngram[0:len(ngram)-1]]
        #print denominator_part
        if denominator_part == 0:
            prob_tilde_list[ngram]=0
        else:
            prob_tilde_list[ngram]=counter_part/(0.0+denominator_part)

    return prob_tilde_list

def alpha(prob_tilde_list,prob_plus1tilde_list,ngram_letter_counter,n_plus1_gram_letter_counter):
    sum_n_1gram_list = 0
    sum_ngram_list = 0

    for ngram in n_plus1_gram_letter_counter:
        if ngram_letter_counter[ngram[0:len(ngram)-1]] > 0:
            sum_ngram_list += prob_plus1tilde_list[ngram]
            sum_n_1gram_list += prob_tilde_list[ngram[1:]]
            #print ngram                
            #print prob_plus1tilde_list[ngram]
            #print ngram[1:]
            #print prob_tilde_list[ngram[1:]]
    #print sum_ngram_list
    #print sum_n_1gram_list
    
    alpha = (1.0-sum_ngram_list)/(1-sum_n_1gram_list)

    return alpha

def bigram_backoff(bigram_letter_counter,unigram_letter_counter,prob_tilde,alpha_n_1gram,all_letter_nr):
    backoff_prob_list = collections.defaultdict(int)

    for ngram,value in bigram_letter_counter.items():
        #print ngram
        if value > 0:
            backoff_prob_list[ngram]=prob_tilde[ngram]
            #print prob_tilde[ngram]
            #print '~'*20
        elif unigram_letter_counter[ngram[1:]] > 0:
            #print alpha_n_1gram
            #print ngram[1:]
            #print unigram_letter_counter[ngram[1:]]
            #print '+'*20
            backoff_prob_list[ngram]=alpha_n_1gram * unigram_letter_counter[ngram[1:]] / (0.0+all_letter_nr)
        
    return backoff_prob_list

def trigram_backoff(trigram_letter_counter,bigram_letter_counter,unigram_letter_counter,prob_tilde,prob_1tilde,bigram_alpha,unigram_alpha,all_leters_nr):
    trigram_backoff=collections.defaultdict(int)
    
    for ngram, value in trigram_letter_counter.items():
        #print ngram
        #print value
        if value > 0:
        #    print 'trigram "'+str(ngram)+'": '+str(value)
        #    print 'prob tilde: '+str(prob_tilde[ngram])
        #    print '---------------------'
            trigram_backoff[ngram]=prob_tilde[ngram]
        elif bigram_letter_counter[ngram[1:]] >0 :
        #    print 'trigram "'+str(ngram)+'": '+str(value)
        #    print 'bigram "'+str(ngram[1:])+'": '+str(bigram_letter_counter[ngram[1:]])
        #    print 'bigram alpha: '+str(bigram_alpha)
        #    print 'unigram prob tilde: '+str(prob_1tilde[ngram[1:]])
        #    print '*********************'
            trigram_backoff[ngram]=bigram_alpha*(0.0+prob_1tilde[ngram[1:]])
        else:
        #    print 'trigram "'+str(ngram)+'": '+str(value)
        #    print 'bigram "'+str(ngram[1:])+'": '+str(bigram_letter_counter[ngram[1:]])
        #    print 'unigram "'+str(ngram[2:])+'": '+str(unigram_letter_counter[ngram[2:]])
        #    print 'unigram aplpha: '+str(unigram_alpha)
        #    print unigram_letter_counter[ngram[2:]]
        #    print tokens_nr
        #    print '++++++++++++++++++++++'
            trigram_backoff[ngram]=unigram_alpha*(0.0+ unigram_letter_counter[ngram[2:]]/tokens_nr)
        
    return trigram_backoff

def four_backoff(four_gram_letter_counter,prob_tilde,trigram_backoff,trigram_alpha):
    four_backoff=collections.defaultdict(int)
    
    for ngram,value in four_gram_letter_counter.items():
        if value > 0:
            four_backoff[ngram]=prob_tilde[ngram]
        else:
            four_backoff[ngram]=(0.0+trigram_alpha)*trigram_backoff[ngram[1:]]
        
    return four_backoff

def fifth_backoff(fifthgram_letter_counter,prob_tilde,four_backoff,fourgram_alpha):
    fifth_backoff=collections.defaultdict(int)
    
    for ngram,value in fifthgram_letter_counter.items():
        if value > 0:
            fifth_backoff[ngram]=prob_tilde[ngram]
        else:
            fifth_backoff[ngram]=(0.0+fourgram_alpha)*four_backoff[ngram[1:]]
        
    return fifth_backoff

#################### REPLACE RARLY OCCURED CHARACTERS WITH * ##################
def list_of_small_nr_of_special_char(unigram_letter_counter):
    list_char=''
    
    for letter in unigram_letter_counter:
        if unigram_letter_counter[letter] <= 30:
            list_char=list_char+letter
    return list_char
def replace_special_char_with_star(filename_in,filenam_out,special_chars):
    data=file(filename_in).read()
    file_out=open(filenam_out,'w')
    i=0
    new_data=data
    
    while i<len(special_chars):
        new_data=new_data.replace(special_chars[i],'*')
        i += 1
    file_out.write(new_data)
    
def tokens_words_with_stars(tokens_list,special_char_string):
    new_tokens_list=collections.defaultdict(int)
    #print special_char_string
    for ngram,value in tokens_list.iteritems():
        new_ngram=ngram
        for c in special_char_string:
            new_ngram=new_ngram.replace(c,'*')
        new_tokens_list[new_ngram]=value
        #print ngram
        #print new_ngram
    return new_tokens_list

        
#################### TEST #######################################################

def test_part(test_tokens_list,order,ngram_prob):
    test_list=defaultdict(Counter)
    result = 0

    for word in test_tokens_list:
        word_result=0
        #print word
        for i in range(len(word)-order+1):
            if ngram_prob[word[i:i+order]] != 0:
            #    print word[i:i+order]
            #    print ngram_prob[word[i:i+order]]
                word_result += math.log(ngram_prob[word[i:i+order]])
        result += word_result
        #print word_result
    return result

def ngram_words(lst,n):
    tlst=lst
    while True:
        a,b = tee(tlst)
        l=tuple(islice(a,n))
        if len(l) == n:
            yield l
            next(b)
            tlst = b
        else:
            break
        
def convert_ngram(ngram_counter):
    ngram_counter_list=collections.defaultdict(int)
    for ngram,value in ngram_counter.iteritems():
        alma=''
        korte=''
        for n in ngram:
            alma=n+ ' '
            korte+=alma
        length=len(korte)
        value2=korte[:len(korte)-1]
        ngram_counter_list[value2]=value
    
    return ngram_counter_list


############################### MAIN ############################
if __name__ == '__main__':
    # LANGUAGE 1
    train_file_name_in1 = '/home/agotaillyes/text_corpus/english.txt'
    language1='english'
    order=2
    
    train_token_list_first1=get_tokens_list(order-1,train_file_name_in1)
    unigram_letter_counter1 = uni_gram_letter_counter(1,train_token_list_first1)
    special_char_list1 = list_of_small_nr_of_special_char(unigram_letter_counter1)    
    train_token_list1=tokens_words_with_stars(train_token_list_first1,special_char_list1)
    
    unigram_letter_counter1 = uni_gram_letter_counter(1,train_token_list1)
    bigram_letter_counter1 = ngram_letter_counter(2,unigram_letter_counter1,train_token_list1)
##    trigram_letter_counter1 = ngram_letter_counter(3,unigram_letter_counter1,train_token_list1)
##    fourgram_letter_counter1 = ngram_letter_counter(4,unigram_letter_counter1,train_token_list1)
##        
    letter_types_nr1 = len(ngram_types(unigram_letter_counter1))
    bigram_types_nr1 = len(ngram_types(bigram_letter_counter1))
##    trigram_types_nr1 = len(ngram_types(trigram_letter_counter1))
##    fourgram_types_nr1 = len(ngram_types(fourgram_letter_counter1)) 
    
    # all letters number
    tokens_nr1 = sum(unigram_letter_counter1.values())
    
    # LANGUAGE 2
    train_file_name_in2 = '/home/agotaillyes/text_corpus/dutch.txt'
    language2='dutch'
    
    train_token_list_first2=get_tokens_list(order-1,train_file_name_in2)
    unigram_letter_counter2 = uni_gram_letter_counter(1,train_token_list_first2)
    special_char_list2 = list_of_small_nr_of_special_char(unigram_letter_counter2)    
    train_token_list2=tokens_words_with_stars(train_token_list_first2,special_char_list2)
    
    unigram_letter_counter2 = uni_gram_letter_counter(1,train_token_list2)
    bigram_letter_counter2 = ngram_letter_counter(2,unigram_letter_counter2,train_token_list2)
##    trigram_letter_counter2 = ngram_letter_counter(3,unigram_letter_counter2,train_token_list2)
##    fourgram_letter_counter2 = ngram_letter_counter(4,unigram_letter_counter2,train_token_list2)
        
    letter_types_nr2 = len(ngram_types(unigram_letter_counter2))
    bigram_types_nr2 = len(ngram_types(bigram_letter_counter2))
##    trigram_types_nr2 = len(ngram_types(trigram_letter_counter2))
##    fourgram_types_nr2 = len(ngram_types(fourgram_letter_counter2)) 
    
    # all letters number
    tokens_nr2 = sum(unigram_letter_counter2.values())
    
    # LANGUAGE 3
    train_file_name_in3 = '/home/agotaillyes/text_corpus/german.txt'
    language3='german'
    
    train_token_list_first3=get_tokens_list(order-1,train_file_name_in3)
    unigram_letter_counter3 = uni_gram_letter_counter(1,train_token_list_first3)
    special_char_list3 = list_of_small_nr_of_special_char(unigram_letter_counter3)    
    train_token_list3=tokens_words_with_stars(train_token_list_first3,special_char_list3)
    
    unigram_letter_counter3 = uni_gram_letter_counter(1,train_token_list3)
    bigram_letter_counter3 = ngram_letter_counter(2,unigram_letter_counter3,train_token_list3)
##    trigram_letter_counter3 = ngram_letter_counter(3,unigram_letter_counter3,train_token_list3)
##    fourgram_letter_counter3 = ngram_letter_counter(4,unigram_letter_counter3,train_token_list3)
        
    letter_types_nr3 = len(ngram_types(unigram_letter_counter3))
    bigram_types_nr3 = len(ngram_types(bigram_letter_counter3))
##    trigram_types_nr3 = len(ngram_types(trigram_letter_counter3))
##    fourgram_types_nr3 = len(ngram_types(fourgram_letter_counter3)) 
    
    # all letters number
    tokens_nr3 = sum(unigram_letter_counter3.values())
    
    #LANGUAGE 4
    train_file_name_in4 = '/home/agotaillyes/text_corpus/danish.txt'
    language4='danish'
    
    train_token_list_first4=get_tokens_list(order-1,train_file_name_in4)
    unigram_letter_counter4 = uni_gram_letter_counter(1,train_token_list_first4)
    special_char_list4 = list_of_small_nr_of_special_char(unigram_letter_counter4)    
    train_token_list4=tokens_words_with_stars(train_token_list_first4,special_char_list4)
    
    unigram_letter_counter4 = uni_gram_letter_counter(1,train_token_list4)
    bigram_letter_counter4 = ngram_letter_counter(2,unigram_letter_counter4,train_token_list4)
##    trigram_letter_counter4 = ngram_letter_counter(3,unigram_letter_counter4,train_token_list4)
##    fourgram_letter_counter4 = ngram_letter_counter(4,unigram_letter_counter4,train_token_list4)
        
    letter_types_nr4 = len(ngram_types(unigram_letter_counter4))
    bigram_types_nr4 = len(ngram_types(bigram_letter_counter4))
##    trigram_types_nr4 = len(ngram_types(trigram_letter_counter4))
##    fourgram_types_nr4 = len(ngram_types(fourgram_letter_counter4)) 
    
    # all letters number
    tokens_nr4 = sum(unigram_letter_counter4.values())
    
##    #LANGUAGE 5
##    train_file_name_in5 = sys.argv[9]
##    language5=sys.argv[10]
##    
##    train_token_list_first5=get_tokens_list(order-1,train_file_name_in5)
##    unigram_letter_counter5 = uni_gram_letter_counter(1,train_token_list_first5)
##    special_char_list5 = list_of_small_nr_of_special_char(unigram_letter_counter5)    
##    train_token_list5=tokens_words_with_stars(train_token_list_first5,special_char_list5)
##    
##    unigram_letter_counter5 = uni_gram_letter_counter(1,train_token_list5)
##    bigram_letter_counter5 = ngram_letter_counter(2,unigram_letter_counter5,train_token_list5)
##    trigram_letter_counter5 = ngram_letter_counter(3,unigram_letter_counter5,train_token_list5)
##    fourgram_letter_counter5 = ngram_letter_counter(4,unigram_letter_counter5,train_token_list5)
##        
##    letter_types_nr5 = len(ngram_types(unigram_letter_counter5))
##    bigram_types_nr5 = len(ngram_types(bigram_letter_counter5))
##    trigram_types_nr5 = len(ngram_types(trigram_letter_counter5))
##    fourgram_types_nr5 = len(ngram_types(fourgram_letter_counter5)) 
##    
##    # all letters number
##    tokens_nr5 = sum(unigram_letter_counter5.values())
    
    test_file_name_in=sys.argv[1]
##    test_token_list = get_tokens_list(order-1,test_file_name_in)
    test_result=collections.defaultdict(int)
    
##    bigram_unsmoothing_prob1 = ngram_unsmoothing_prob(2,bigram_letter_counter1,unigram_letter_counter1,tokens_nr1)
##    bigram_unsmoothing_prob2 = ngram_unsmoothing_prob(2,bigram_letter_counter2,unigram_letter_counter2,tokens_nr2)
##    bigram_unsmoothing_prob3 = ngram_unsmoothing_prob(2,bigram_letter_counter3,unigram_letter_counter3,tokens_nr3)
##    bigram_unsmoothing_prob4 = ngram_unsmoothing_prob(2,bigram_letter_counter4,unigram_letter_counter4,tokens_nr4)
##    bigram_unsmoothing_prob5 = ngram_unsmoothing_prob(2,bigram_letter_counter5,unigram_letter_counter5,tokens_nr5)
    
##    bigram_add_one_prob1 = ngram_add_one_prob(2,bigram_letter_counter1,unigram_letter_counter1,tokens_nr1,letter_types_nr1)
##    bigram_add_one_prob2 = ngram_add_one_prob(2,bigram_letter_counter2,unigram_letter_counter2,tokens_nr2,letter_types_nr2)
##    bigram_add_one_prob3 = ngram_add_one_prob(2,bigram_letter_counter3,unigram_letter_counter3,tokens_nr3,letter_types_nr3)
##    bigram_add_one_prob4 = ngram_add_one_prob(2,bigram_letter_counter4,unigram_letter_counter4,tokens_nr4,letter_types_nr4)
##    bigram_add_one_prob5 = ngram_add_one_prob(2,bigram_letter_counter5,unigram_letter_counter5,tokens_nr5,letter_types_nr5)
    
##    bigram_witten_bell1 = ngram_witten_bell_prob(2,bigram_letter_counter1,tokens_nr1,letter_types_nr1,letter_types_nr1)
##    bigram_witten_bell2 = ngram_witten_bell_prob(2,bigram_letter_counter2,tokens_nr2,letter_types_nr2,letter_types_nr2)
##    bigram_witten_bell3 = ngram_witten_bell_prob(2,bigram_letter_counter3,tokens_nr3,letter_types_nr3,letter_types_nr3)
##    bigram_witten_bell4 = ngram_witten_bell_prob(2,bigram_letter_counter4,tokens_nr4,letter_types_nr4,letter_types_nr4)
##    bigram_witten_bell5 = ngram_witten_bell_prob(2,bigram_letter_counter5,tokens_nr5,letter_types_nr5,letter_types_nr5)
##    
##    bigram_star1 = c_star(bigram_letter_counter1,bigram_witten_bell1,tokens_nr1)
##    bigram_star2 = c_star(bigram_letter_counter2,bigram_witten_bell2,tokens_nr2)
##    bigram_star3 = c_star(bigram_letter_counter3,bigram_witten_bell3,tokens_nr3)
##    bigram_star4 = c_star(bigram_letter_counter4,bigram_witten_bell4,tokens_nr4)
##    bigram_star5 = c_star(bigram_letter_counter5,bigram_witten_bell5,tokens_nr5)
##    
##    unigram_prob_tilde1 = unigram_prob_tilde(unigram_letter_counter1,tokens_nr1)
##    unigram_prob_tilde2 = unigram_prob_tilde(unigram_letter_counter2,tokens_nr2)
##    unigram_prob_tilde3 = unigram_prob_tilde(unigram_letter_counter3,tokens_nr3)
##    unigram_prob_tilde4 = unigram_prob_tilde(unigram_letter_counter4,tokens_nr4)
##    unigram_prob_tilde5 = unigram_prob_tilde(unigram_letter_counter5,tokens_nr5)
##    
##    bigram_prob_tilde1 = prob_tilde(bigram_letter_counter1,unigram_letter_counter1,bigram_star1)
##    bigram_prob_tilde2 = prob_tilde(bigram_letter_counter2,unigram_letter_counter2,bigram_star2)
##    bigram_prob_tilde3 = prob_tilde(bigram_letter_counter3,unigram_letter_counter3,bigram_star3)
##    bigram_prob_tilde4 = prob_tilde(bigram_letter_counter4,unigram_letter_counter4,bigram_star4)
##    bigram_prob_tilde5 = prob_tilde(bigram_letter_counter5,unigram_letter_counter5,bigram_star5)
##
##    unigram_alpha1 = alpha(unigram_prob_tilde1,bigram_prob_tilde1,unigram_letter_counter1,bigram_letter_counter1)
##    unigram_alpha2 = alpha(unigram_prob_tilde2,bigram_prob_tilde2,unigram_letter_counter2,bigram_letter_counter2)
##    unigram_alpha3 = alpha(unigram_prob_tilde3,bigram_prob_tilde3,unigram_letter_counter3,bigram_letter_counter3)
##    unigram_alpha4 = alpha(unigram_prob_tilde4,bigram_prob_tilde4,unigram_letter_counter4,bigram_letter_counter4)
##    unigram_alpha5 = alpha(unigram_prob_tilde5,bigram_prob_tilde5,unigram_letter_counter5,bigram_letter_counter5)
##    
##    bigram_backoff_prob1=bigram_backoff(bigram_letter_counter1,unigram_letter_counter1,bigram_prob_tilde1,unigram_alpha1,tokens_nr1)
##    bigram_backoff_prob2=bigram_backoff(bigram_letter_counter2,unigram_letter_counter2,bigram_prob_tilde2,unigram_alpha2,tokens_nr2)
##    bigram_backoff_prob3=bigram_backoff(bigram_letter_counter3,unigram_letter_counter3,bigram_prob_tilde3,unigram_alpha3,tokens_nr3)
##    bigram_backoff_prob4=bigram_backoff(bigram_letter_counter4,unigram_letter_counter4,bigram_prob_tilde4,unigram_alpha4,tokens_nr4)
##    bigram_backoff_prob5=bigram_backoff(bigram_letter_counter5,unigram_letter_counter5,bigram_prob_tilde5,unigram_alpha5,tokens_nr5)

    bigram_add_one_prob1 = ngram_add_one_prob(2,bigram_letter_counter1,unigram_letter_counter1,tokens_nr1,letter_types_nr1)
    bigram_add_one_prob2 = ngram_add_one_prob(2,bigram_letter_counter2,unigram_letter_counter2,tokens_nr2,letter_types_nr2)
    bigram_add_one_prob3 = ngram_add_one_prob(2,bigram_letter_counter3,unigram_letter_counter3,tokens_nr3,letter_types_nr3)
    bigram_add_one_prob4 = ngram_add_one_prob(2,bigram_letter_counter4,unigram_letter_counter4,tokens_nr4,letter_types_nr4)

###################################### TEST PART ###########################################
    with open(test_file_name_in,"r") as myfile:
        data=myfile.read()
    words_test = re.findall("\w+",data)
    
    bigram_counter_test = Counter(ngram_words(words_test,2))
    bigram_letter_counter_test =convert_ngram(bigram_counter_test)
    
    counter=collections.defaultdict(int)
    
    for ngram,value in bigram_letter_counter_test.iteritems():
        actual_list=collections.defaultdict(int)
        for word in ngram.split():
            actual_list[word]+=1
        result1=test_part(actual_list,order,bigram_add_one_prob1)
        result2=test_part(actual_list,order,bigram_add_one_prob2)
        result3=test_part(actual_list,order,bigram_add_one_prob3)
        result4=test_part(actual_list,order,bigram_add_one_prob4)
        
        test_result[language1]=result1
        test_result[language2]=result2
        test_result[language3]=result3
        test_result[language4]=result4
        
        sorted_test = sorted(test_result.iteritems(), key=lambda (k,v):v,reverse=True)
        counter[sorted_test[0][0]]+=1
    print counter