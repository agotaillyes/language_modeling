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
def unigram_letter_counter(j,tokens_list):
    ngram_list = collections.defaultdict(int)

    for word in tokens_list:
        lenght = len(word)
        i = 0
        while lenght-j+1 > i:
            ngram_list[word[i:i+j]] += tokens_list[word]
            i += 1
    space_value=ngram_list[' ']
    space_value=(space_value-2*(j-1))/2
    ngram_list[' ']=space_value

    return ngram_list

def ngram_letter_counter(j,unigram_letter_counter,tokens_list):
    vocabulary=[]
    ngram_letter_counter = collections.defaultdict(int)
    
    for unigram,value in unigram_letter_counter.iteritems():
        vocabulary.append(unigram);
    
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
    
    for ngram,value in ngram_letter_counter.iteritems():
        if n_1gram_letter_counter[ngram[0:i-1]] != 0:
            ngram_prob_list[ngram] = (value+0.0) / n_1gram_letter_counter[ngram[0:i-1]]
        else:
            ngram_prob_list[ngram] = 0.0
        
    return ngram_prob_list

#################### REPLACE RARLY OCCURED CHARACTERS WITH * ##################
def list_of_small_nr_of_special_char(unigram_letter_counter):
    list_char=''
    
    for letter,value in unigram_letter_counter.iteritems():
        if value <= 30:
            list_char=list_char+letter
    return list_char
    
def tokens_words_with_stars(tokens_list,special_char_string):
    new_tokens_list=collections.defaultdict(int)
    for ngram,value in tokens_list.iteritems():
        new_ngram=ngram
        for c in special_char_string:
            new_ngram=new_ngram.replace(c,'*')
        new_tokens_list[new_ngram]=value
    return new_tokens_list

#################### TEST #######################################################

def test_part(test_tokens_list,order,ngram_prob):
    test_list=defaultdict(Counter)
    result = 0

    for word in test_tokens_list:
        word_result=0
        for i in range(len(word)-order+1):
            if ngram_prob[word[i:i+order]] != 0:
                word_result += math.log(ngram_prob[word[i:i+order]])
        result += word_result
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
            if not n.startswith('-') and n != '':
                n = re.sub('[!"#$%&\'()+,-./:;<=>?@[\\]^_`{|}~1234567890]','',n)
                n = n.lower()
                alma=n+ ' '
                korte+=alma
        length=len(korte)
        value2=korte[:len(korte)-1]
        ngram_counter_list[value2]=value
    
    return ngram_counter_list

    ngram_counter_list=collections.defaultdict(int)
    for ngram,value in ngram_counter.iteritems():
        alma=''
        korte=''
        for n in ngram:
            if not n.startswith('-') and n != '':
                n = re.sub('[!"#$%&\'()+,-./:;<=>?@[\\]^_`{|}~1234567890]','',n)
                n = n.lower()
                alma=n+ ' '
                korte+=alma
        length=len(korte)
        value2=korte[:len(korte)-1]
        ngram_counter_list[value2]=value
    
    return ngram_counter_list

############################### MAIN ############################
if __name__ == '__main__':
    # LANGUAGE 1
    train_file_name_in1 = '/home/agotaillyes/text_corpus/english2.txt'
    language1='english'
    order=2
    
    train_token_list_first1=get_tokens_list(order-1,train_file_name_in1)
    unigram_letter_counter1 = unigram_letter_counter(1,train_token_list_first1)
    special_char_list1 = list_of_small_nr_of_special_char(unigram_letter_counter1)    
    train_token_list1=tokens_words_with_stars(train_token_list_first1,special_char_list1)
    
    unigram_letter_counter1 = unigram_letter_counter(1,train_token_list1)
    bigram_letter_counter1 = ngram_letter_counter(2,unigram_letter_counter1,train_token_list1)
        
    letter_types_nr1 = len(ngram_types(unigram_letter_counter1))
    
    # all letters number
    tokens_nr1 = sum(unigram_letter_counter1.values())
    
    # LANGUAGE 2
    train_file_name_in2 = '/home/agotaillyes/text_corpus/dutch2.txt'
    language2='dutch'
    
    train_token_list_first2=get_tokens_list(order-1,train_file_name_in2)
    unigram_letter_counter2 = unigram_letter_counter(1,train_token_list_first2)
    special_char_list2 = list_of_small_nr_of_special_char(unigram_letter_counter2)    
    train_token_list2=tokens_words_with_stars(train_token_list_first2,special_char_list2)
    
    unigram_letter_counter2 = unigram_letter_counter(1,train_token_list2)
    bigram_letter_counter2 = ngram_letter_counter(2,unigram_letter_counter2,train_token_list2)
        
    letter_types_nr2 = len(ngram_types(unigram_letter_counter2))
    
    # all letters number
    tokens_nr2 = sum(unigram_letter_counter2.values())
    
    # LANGUAGE 3
    train_file_name_in3 = '/home/agotaillyes/text_corpus/german2.txt'
    language3='german'
    
    train_token_list_first3=get_tokens_list(order-1,train_file_name_in3)
    unigram_letter_counter3 = unigram_letter_counter(1,train_token_list_first3)
    special_char_list3 = list_of_small_nr_of_special_char(unigram_letter_counter3)    
    train_token_list3=tokens_words_with_stars(train_token_list_first3,special_char_list3)
    
    unigram_letter_counter3 = unigram_letter_counter(1,train_token_list3)
    bigram_letter_counter3 = ngram_letter_counter(2,unigram_letter_counter3,train_token_list3)
        
    letter_types_nr3 = len(ngram_types(unigram_letter_counter3))
    
    # all letters number
    tokens_nr3 = sum(unigram_letter_counter3.values())
    
    #LANGUAGE 4
    train_file_name_in4 = '/home/agotaillyes/text_corpus/danish2.txt'
    language4='danish'
    
    train_token_list_first4=get_tokens_list(order-1,train_file_name_in4)
    unigram_letter_counter4 = unigram_letter_counter(1,train_token_list_first4)
    special_char_list4 = list_of_small_nr_of_special_char(unigram_letter_counter4)    
    train_token_list4=tokens_words_with_stars(train_token_list_first4,special_char_list4)
    
    unigram_letter_counter4 = unigram_letter_counter(1,train_token_list4)
    bigram_letter_counter4 = ngram_letter_counter(2,unigram_letter_counter4,train_token_list4)
        
    letter_types_nr4 = len(ngram_types(unigram_letter_counter4))
    bigram_types_nr4 = len(ngram_types(bigram_letter_counter4))
    
    # all letters number
    tokens_nr4 = sum(unigram_letter_counter4.values())
    
    test_file_name_in=sys.argv[1]
    test_result=collections.defaultdict(int)
    
    bigram_unsmoothing_prob1 = ngram_unsmoothing_prob(2,bigram_letter_counter1,unigram_letter_counter1,tokens_nr1)
    bigram_unsmoothing_prob2 = ngram_unsmoothing_prob(2,bigram_letter_counter2,unigram_letter_counter2,tokens_nr2)
    bigram_unsmoothing_prob3 = ngram_unsmoothing_prob(2,bigram_letter_counter3,unigram_letter_counter3,tokens_nr3)
    bigram_unsmoothing_prob4 = ngram_unsmoothing_prob(2,bigram_letter_counter4,unigram_letter_counter4,tokens_nr4)

###################################### TEST PART ###########################################
    with open(test_file_name_in,"r") as myfile:
        data=myfile.read()
    words_test = re.findall("\w+",data)
    ngram_test=sys.argv[2]
    ngram_test_nr=int(ngram_test,0)
##    print 'ngram: ' + str(ngram_test_nr)
    bigram_counter_test = Counter(ngram_words(words_test,ngram_test_nr))
    bigram_letter_counter_test =convert_ngram(bigram_counter_test)
##    print 'bigram test german'
##    print sum(bigram_letter_counter_test.values())
    
    counter=collections.defaultdict(int)
    
    for ngram,value in bigram_letter_counter_test.iteritems():
        actual_list=collections.defaultdict(int)
        for word in ngram.split():
            actual_list[word]+=1
        result1=test_part(actual_list,order,bigram_unsmoothing_prob1)
        result2=test_part(actual_list,order,bigram_unsmoothing_prob2)
        result3=test_part(actual_list,order,bigram_unsmoothing_prob3)
        result4=test_part(actual_list,order,bigram_unsmoothing_prob4)
        
        test_result[language1]=result1
        test_result[language2]=result2
        test_result[language3]=result3
        test_result[language4]=result4
        
        sorted_test = sorted(test_result.iteritems(), key=lambda (k,v):v,reverse=True)
        counter[sorted_test[0][0]]+=1
    
    norm_nr = sum(counter.values())
    for lang, value in counter.iteritems():
        counter[lang]=(value+0.0)/norm_nr
        
    result = sorted(counter.iteritems(), key=lambda (k,v):v,reverse=True)
    print result
    
##    print 'unsmoothing'
##    print sorted(bigram_unsmoothing_prob1.iteritems(), key=lambda (k,v):v,reverse=True)