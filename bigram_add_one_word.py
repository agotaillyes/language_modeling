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
import psutil
import time

# megkapjuk a beolvasott file osszes token-et (szavat) kozpontozas nelkul
def get_tokens_list(input_file_name):
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
                    if word != '':
                        tokens_list[word] += 1
    return tokens_list
 
 # osszeszamolja, hogy j egymast koveto karakterbol mennyi van
 # a szovegben
         
def ngram(lst,n):
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
    szilva=collections.defaultdict(int)
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
        value2 = korte[:len(korte)-1]
        szilva[value2]=value
    return szilva
 
def n_gram_word_counter(j,ngram_word_counter,tokens_list):
    ngram_counter=collections.defaultdict(int)
    if j==1:
        ngram_l = (item for item,value in tokens_list.iteritems())
    else:
        ngram_l = (' '.join(item) for item in itertools.product(tokens_list,repeat=j))
     
    for ngram in ngram_l:
        ngram_counter[ngram]=0
  
    for ngram,values in ngram_word_counter.iteritems():
        ngram_counter[ngram]=values
 
    return ngram_counter
    
#################### ADD-ONE probabilities ############################

def ngram_add_one_prob(i,ngram_letter_counter,n_1gram_letter_counter,tokens_nr,letter_types_nr):
    ngram_prob_list = collections.defaultdict(int)
    
    for ngram,value in ngram_letter_counter.iteritems():
        list=ngram.split()
        first_word=list[0]
        ngram_prob_list[ngram] = (value+1.0) / (n_1gram_letter_counter[first_word]+letter_types_nr)
        
    return ngram_prob_list

########################### TRAIN #############################
def train_word_ngram(nplus1_gram_word_prob):
    lm = defaultdict(Counter)
 
    for word,value in nplus1_gram_word_prob.iteritems():
        if value!=0:
##            print word
            data = word
            list=word.split()
            history=list[0]
            char=list[1]
##            print value
##            print history
##            print char
##            print '~'*80               
            lm[history][char]=value
##    print lm
    return lm
 
def convert_train_to_rank(train_ngram_prob):
    train_rank=defaultdict(Counter)
    
    for keys, values in train_ngram_prob.iteritems():
        train_list=[keys1 for values1,keys1 in values.items()]
        rank_list = rankdata(train_list,method='max')
        reversed_rank_list=[]
        for element in rank_list:
            alma = len(rank_list)+1-element
            reversed_rank_list.append(alma)
        for index,item in enumerate(values):
            train_rank[keys][item]=reversed_rank_list[index]
    return train_rank
 
def list_of_small_nr_of_word(unigram_word_counter_train):

    list_word = set()
     
    for word,value in unigram_word_counter_train.iteritems():
         if value == 1:
            list_word.add(word)
    return list_word

def replace_special_word_with_rare(filename_in,filenam_out,rare_words):
    with open(filename_in,'r') as infile, open(filenam_out,'w') as outfile:
        for line in infile:
            for word in line.split():
                if not word.startswith('-') and word != '':
                    word=word.lower()
                    word = re.sub('[!"#$%&\'()+,-./:;<=>?@[\\]^_`{|}~1234567890]','',word)
                    if word.lower() in rare_words:
                        word = '_RARE_'
                    outfile.write(word+' ')
        
#################### TEST #######################################################

def test_part(ngram_letter_counter):
    test_list=defaultdict(Counter)
 
    for word,value in ngram_letter_counter.iteritems():
        if value !=0:
            list = word.split()
            history=list[0]
            char=list[1]
            alma=word
            test_list[history][char] = value
##            print alma
##            print test_list[history][char]
##            print '~'*80
##    print test_list
    return test_list

def list_of_oov_words(train_token_list,test_token_list):
    list_words=set()
     
    for ngram in test_token_list:
        if ngram not in train_token_list:
            list_words.add(ngram)
             
    return list_words
        
################################ PRINT RESULTS ################################
def result_list(test_counter, ranked_train_counter,normalize_nr):
    result =defaultdict(Counter)
    normalized_result =defaultdict(Counter)
     
    i=1
    while i < 11:
        result[i]=0
        i += 1
     
    for ngram,test_list in test_counter.iteritems():
##        print '~' * 80 
##        print 'test ngram: ' + str(ngram) + ' - list: '+str(test_list)
        list=test_list
        if len(ranked_train_counter[ngram]) > 10:
            train_list=sorted(ranked_train_counter[ngram].iteritems(), key=lambda x:x[1])[:10]
        else:
            train_list=sorted(ranked_train_counter[ngram].iteritems(), key=lambda x:x[1])
        array=[k for k,v in train_list]
##        print 'train list: ' + str(train_list)
        for key,value in list.items():
##            print 'list element: ' +str(array) + ' - ' +str(key)+' : '+str(value)
            if key in array:
                korte=value
                for train_key,train_values in train_list:
                    if key == train_key:
##                        print 'train key-test key: '+str(train_key) +'-'+str(key)
                        alma=int(float(train_values))
##                        print 'train value: '+str(alma)
                        for i in xrange(alma,len(result)+1,1):
                            result[i]+=korte
##                            print 'result['+str(i)+']='+str(result[i])
##                        print '~' *80
     
    for key,value in result.items():
        value_new =(value+0.0) / normalize_nr
        normalized_result[key]=value_new
    return normalized_result

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem
         
start_time = time.time()
############################### MAIN ############################
if __name__ == '__main__':
    train_file_name_in = sys.argv[1]
    train_file_name_out = 'train.txt'
    test_file_name_in=sys.argv[2]
    test_file_name_out='test.txt'
         
    train_token_list_first=get_tokens_list(train_file_name_in)
    test_token_list_first=get_tokens_list(test_file_name_in)
     
    rare_words = list_of_small_nr_of_word(train_token_list_first)
    replace_special_word_with_rare(train_file_name_in,train_file_name_out,rare_words)
    train_token_list = get_tokens_list(train_file_name_out)
     
    list_of_oov_in_test_file = list_of_oov_words(train_token_list,test_token_list_first)
    replace_special_word_with_rare(test_file_name_in,test_file_name_out,list_of_oov_in_test_file)
    test_token_list=get_tokens_list(test_file_name_out)
    test_tokens_nr = sum(test_token_list.values())
     
    tokens_nr = sum(train_token_list.values())
    unigram_word_counter_train = train_token_list
##    print sorted(unigram_word_counter_train.iteritems(),key=lambda (k,v):v,reverse=True)
    word_type = len(unigram_word_counter_train)
    #print sorted(test_token_list.iteritems(),key=lambda (k,v):v,reverse=True)
     
    with open(train_file_name_out,"r") as myfile:
        data=myfile.read()
    words_train = re.findall("\w+",data)
     
    bigram_counter_train = Counter(ngram(words_train,2))
    bigram_count_train = convert_ngram(bigram_counter_train)
    bigram_word_counter_train=n_gram_word_counter(2,bigram_count_train,train_token_list)
##    print sorted(unigram_word_counter_train.iteritems(),key=lambda (k,v):v,reverse=True)
##    print 'token nr: ' + str(tokens_nr)
##    print 'word type: ' + str(word_type)
     
    with open(test_file_name_out,"r") as myfile:
        data=myfile.read()
    words_test = re.findall("\w+",data)
     
    bigram_counter_test = Counter(ngram(words_test,2))
    bigram_count_test = convert_ngram(bigram_counter_test)
    bigram_word_counter_test=n_gram_word_counter(2,bigram_count_test,test_token_list)
     
    bigram_add_one_prob = ngram_add_one_prob(2,bigram_word_counter_train,unigram_word_counter_train,tokens_nr,word_type) 
    
    train_unigram_add_one = train_word_ngram(bigram_add_one_prob)    
    test_unigram_add_one = test_part(bigram_word_counter_test) 
    
    ranked_train = convert_train_to_rank(train_unigram_add_one)
    
    all_test_bigrams_nr = sum(bigram_word_counter_test.values())
    
    result= result_list(test_unigram_add_one,ranked_train,all_test_bigrams_nr)
    
##    print 'add-one smoothing'
##    print sorted(bigram_add_one_prob.iteritems(),key=lambda (k,v):v,reverse=True)
##    print sorted(test_unigram_add_one.iteritems(),key=lambda (k,v):v,reverse=False)
 
##    print '~' * 80
##    print 'bigram add-one probability'
##    print 'train word number: ' + str(tokens_nr)
##    print 'test word number: ' + str(test_tokens_nr)
##    for key,value in result.items():
##        print 'top' +str(key)+': ' + str(value)+'%'

    memory=memory_usage_psutil()
    print memory*1.04858
    print("--- %s seconds ---" % (time.time() - start_time))