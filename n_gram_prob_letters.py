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
def uni_gram_letter_counter(j,tokens_list):
    ngram_list = collections.defaultdict(int)

    for word in tokens_list:
        lenght = len(word)
        i = 0
        while lenght-j+1 > i:
            ngram_list[word[i:i+j]] += tokens_list[word]
            i += 1

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

def ngram_witten_bell_prob(i,ngram_letter_counter,tokens_nr,ngram_type_nr,letter_types):
    z_list = collections.defaultdict(int)
    observed_type = collections.defaultdict(int)
    ngram_prob_list = collections.defaultdict(int)
    z = 0

    for ngram1, value1 in ngram_letter_counter.iteritems():
        if i!= 1:
            alma=0
        
            for ngram2,value2 in ngram_letter_counter.iteritems():
                if ngram1[0:i-1]==ngram2[0:i-1] and value2 !=0: 
                    #print str(ngram1)+' '+str(value1)
                    #print str(ngram2)+' '+str(value2)
                    #print ngram1[0:i-1]
                    #print ngram2[0:i-1]
                    alma += 1
           
            observed_type[ngram1]=alma
            z_list[ngram1]=letter_types-alma
            #print alma
            #print letter_types-alma
            #print '~'*80
        else:
            if(value1 == 0):
                z += 1

    for ngram,value in ngram_letter_counter.iteritems():
        if i != 1:
            if(value == 0):
                ngram_prob_list[ngram] = (observed_type[ngram]+0.0)/(z_list[ngram]*(tokens_nr+observed_type[ngram]))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0)/(tokens_nr+observed_type[ngram])
        else:
            if(ngram_letter_counter[ngram] == 0):
                ngram_prob_list[ngram] = (ngram_type_nr + 0.0) / (z * (tokens_nr + ngram_type_nr))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram] + 0.0) / (tokens_nr + ngram_type_nr)
        
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

########################### TRAIN #############################
def train_char_ngram(ngram_letter_counter,order,orderplus1_ngram_prob):
    lm = defaultdict(Counter)

    for word in ngram_letter_counter:
        #print word
        data = word
        #print '~'*80
        history,char =data[0:order],data[order]
        #print history
        #print char
        if char != " " and not (history.startswith(" ") or history.endswith(" ")):                
            word=history+char
            lm[history][char]=orderplus1_ngram_prob[word]
            #print word
            #print lm[history][char]
        #print '~'*80
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
        
#################### TEST #######################################################

def test_part(ngram_letter_counter,order):
    test_list=defaultdict(Counter)

    for word in ngram_letter_counter:
        data=word
        #print '~'*80
        history, char= data[0:order],data[order]
        #print history
        #print char
        if char != " " and not (history.startswith(" ") or history.endswith(" ")):
            alma=history+char
            test_list[history][char] = ngram_letter_counter[alma]
            #print alma
            #print test_list[history][char]
        #print '~'*80
    return test_list

################################ PRINT RESULTS ################################
def result_list(test_counter, ranked_train_counter,normalize_nr):
    result =defaultdict(Counter)
    normalized_result =defaultdict(Counter)
    
    i=1
    while i < 11:
        result[i]=0
        i += 1
    
    for ngram,test_list in test_counter.iteritems():
        #print '~' * 80 
        #print 'test ngram: ' + str(ngram) + ' - list: '+str(test_list)
        list=test_list
        if len(ranked_train_counter[ngram]) > 10:
            train_list=sorted(ranked_train_counter[ngram].iteritems(), key=lambda x:x[1])[:10]
        else:
            train_list=sorted(ranked_train_counter[ngram].iteritems(), key=lambda x:x[1])
        array=[k for k,v in train_list]
        #print 'train list: ' + str(train_list)
        for key,value in list.items():
            #print 'array: ' +str(array) + ' - ' +str(key)+' : '+str(value)
            if key in array:
                korte=value
                for train_key,train_values in train_list:
                    if key == train_key:
                        #print 'train key-test key: '+str(train_key) +'-'+str(key)
                        alma=int(float(train_values))
                        #print 'train value: '+str(alma)
                        for i in xrange(alma,len(result)+1,1):
                            result[i]+=korte
                            #print 'result['+str(i)+']='+str(result[i])
                        #print '~' *80
    
    for key,value in result.items():
        value_new =(value+0.0) / normalize_nr
        normalized_result[key]=value_new
    return normalized_result
        

############################### MAIN ############################
if __name__ == '__main__':
    train_file_name_in = sys.argv[1]
    test_file_name_in=sys.argv[2]
    test_file_name_out=sys.argv[3]
    language=sys.argv[4]
    
    train_token_list_first=get_tokens_list(train_file_name_in)
    unigram_letter_counter = uni_gram_letter_counter(1,train_token_list_first)
    special_char_list = list_of_small_nr_of_special_char(unigram_letter_counter)    
    train_token_list=tokens_words_with_stars(train_token_list_first,special_char_list)

    test_token_list = get_tokens_list(test_file_name_out)


#############################################################################################
    
    unigram_letter_counter = uni_gram_letter_counter(1,train_token_list)
    #all letters number
    tokens_nr = sum(unigram_letter_counter.values())

    bigram_letter_counter = ngram_letter_counter(2,unigram_letter_counter,train_token_list)
    trigram_letter_counter = ngram_letter_counter(3,unigram_letter_counter,train_token_list)
    fourgram_letter_counter = ngram_letter_counter(4,unigram_letter_counter,train_token_list)
    #fifthgram_letter_counter = ngram_letter_counter(5,unigram_letter_counter,train_token_list)
    
    unigram_types = ngram_types(unigram_letter_counter)
    bigram_types = ngram_types(bigram_letter_counter)
    trigram_types = ngram_types(trigram_letter_counter)
    fourgram_types = ngram_types(fourgram_letter_counter)
    
    letter_types_nr = len(unigram_types)
    bigram_types_nr = len(bigram_types)
    trigram_types_nr = len(trigram_types)
    fourgram_types_nr = len(fourgram_types)    
    
    #print len(unigram_letter_counter)
    #print letter_types_nr
    #print len(bigram_letter_counter)
    #print bigram_types_nr
    #print len(trigram_letter_counter)
    #print trigram_types_nr
    
    #all_train_bigrams_nr=sum(bigram_letter_counter.values())
    #all_train_trigram_nr=sum(trigram_letter_counter.values())
    #all_train_fourgram_nr=sum(fourgram_letter_counter.values())
    
    #unigram_unsmoothing_prob = ngram_unsmoothing_prob(1,unigram_letter_counter,unigram_letter_counter,tokens_nr)
    #unigram_add_one_prob = ngram_add_one_prob(1,unigram_letter_counter,unigram_letter_counter,tokens_nr,letter_types_nr)  
    #unigram_witten_bell = ngram_witten_bell_prob(1,unigram_letter_counter,tokens_nr,letter_types_nr,letter_types_nr)
    bigram_unsmoothing_prob = ngram_unsmoothing_prob(2,bigram_letter_counter,unigram_letter_counter,tokens_nr)
    bigram_add_one_prob = ngram_add_one_prob(2,bigram_letter_counter,unigram_letter_counter,tokens_nr,letter_types_nr)  
    bigram_witten_bell = ngram_witten_bell_prob(2,bigram_letter_counter,tokens_nr,letter_types_nr,letter_types_nr)
    #trigram_unsmoothing = ngram_unsmoothing_prob(3,trigram_letter_counter,bigram_letter_counter,tokens_nr)
    #trigram_add_one_prob = ngram_add_one_prob(3,trigram_letter_counter,bigram_letter_counter,tokens_nr,letter_types_nr)
    #trigram_witten_bell = ngram_witten_bell_prob(3,trigram_letter_counter,tokens_nr,trigram_types_nr,letter_types_nr)
    #fourgram_unsmoothing = ngram_unsmoothing_prob(4,fourgram_letter_counter,trigram_letter_counter,tokens_nr)
    #fourgram_add_one_prob = ngram_add_one_prob(4,fourgram_letter_counter,trigram_letter_counter,tokens_nr,letter_types_nr)
    #fourgram_witten_bell = ngram_witten_bell_prob(4,fourgram_letter_counter,tokens_nr,fourgram_types_nr,letter_types_nr)
    #fifth_witten_bell = ngram_witten_bell_prob(5,fifthgram_letter_counter,fourgram_letter_counter)
    #print sorted(trigram_unsmoothing.iteritems(), key=lambda (k,v):v,reverse=False)

    #unigram_star = c_star(unigram_letter_counter,unigram_witten_bell,tokens_nr)
    bigram_star = c_star(bigram_letter_counter,bigram_witten_bell,tokens_nr)
    #trigram_star = c_star(trigram_letter_counter,trigram_witten_bell,tokens_nr)
    #fourgram_star = c_star(fourgram_letter_counter,fourgram_witten_bell,tokens_nr)
    #fifthgram_star = c_star(fifthgram_letter_counter)
   
    unigram_prob_tilde = unigram_prob_tilde(unigram_letter_counter,tokens_nr)
    bigram_prob_tilde = prob_tilde(bigram_letter_counter,unigram_letter_counter,bigram_star)
    #trigram_prob_tilde = prob_tilde(trigram_letter_counter,bigram_letter_counter,trigram_star)
    #fourgram_prob_tilde = prob_tilde(fourgram_letter_counter,trigram_letter_counter,fourgram_star)
    #print sorted(trigram_prob_tilde.iteritems(),key=lambda (k,v):v,reverse=False)
    #fifth_prob_tilde = prob_tilde(fifthgram_letter_counter,fourgram_letter_counter,fifthgram_star)
    
    unigram_alpha = alpha(unigram_prob_tilde,bigram_prob_tilde,unigram_letter_counter,bigram_letter_counter)
    #bigram_alpha = alpha(bigram_prob_tilde,trigram_prob_tilde,bigram_letter_counter,trigram_letter_counter)
    #trigram_alpha = alpha(trigram_prob_tilde,fourgram_prob_tilde,trigram_letter_counter,bigram_letter_counter)
    #fourgram_alpha = alpha(fourgram_prob_tilde,fifth_prob_tilde,fourgram_letter_counter,trigram_letter_counter)

    bigram_backoff_prob=bigram_backoff(bigram_letter_counter,unigram_letter_counter,bigram_prob_tilde,unigram_alpha,tokens_nr)
    #trigram_backoff_prob=trigram_backoff(trigram_letter_counter,bigram_letter_counter,unigram_letter_counter,trigram_prob_tilde,bigram_prob_tilde,bigram_alpha,unigram_alpha,tokens_nr)
    #print sorted(trigram_backoff_prob.iteritems(),key=lambda (k,v):v,reverse=True)
    #four_backoff_prob = four_backoff(fourgram_letter_counter,fourgram_prob_tilde,trigram_backoff_prob,trigram_alpha)
    #fifth_backoff_prob = fifth_backoff(fifthgram_letter_counter,fifth_prob_tilde,four_backoff,fourgram_alpha)
    
    #train_unigram_unsmoothing = train_char_ngram(bigram_letter_counter,1,bigram_unsmoothing_prob)    
    #test_unigram_unsmoothing = test_part(bigram_letter_counter,1)   
    #train_unigram_add_one = train_char_ngram(bigram_letter_counter,1,bigram_add_one_prob)    
    #test_unigram_add_one = test_part(bigram_letter_counter,1) 
    #train_unigram_witten_bell = train_char_ngram(bigram_letter_counter,1,bigram_witten_bell)    
    #test_unigram_witten_bell = test_part(bigram_letter_counter,1)   
    #train_unigram_backoff = train_char_ngram(bigram_letter_counter,1,bigram_backoff_prob)    
    #test_unigram_backoff = test_part(bigram_letter_counter,1) 
    #train_bigram_unsmoothing = train_char_ngram(trigram_letter_counter,2,trigram_unsmoothing)
    #test_bigram_unsmoothing = test_part(trigram_letter_counter,2)
    #train_bigram_add_one = train_char_ngram(trigram_letter_counter,2,trigram_add_one_prob)
    #test_bigram_add_one = test_part(trigram_letter_counter,2)
    #train_bigram_witten_bell = train_char_ngram(trigram_letter_counter,2,trigram_witten_bell)
    #test_bigram_witten_bell = test_part(trigram_letter_counter,2)
    #train_bigram_backoff = train_char_ngram(trigram_letter_counter,2,trigram_backoff_prob)
    #test_bigram_backoff = test_part(trigram_letter_counter,2)
    #train_trigram_unsmoothing = train_char_ngram(fourgram_letter_counter,3,fourgram_unsmoothing)
    #test_trigram_unsmoothing = test_part(fourgram_letter_counter,3)
    #train_trigram_add_one = train_char_ngram(fourgram_letter_counter,3,fourgram_add_one_prob)
    #test_trigram_add_one = test_part(fourgram_letter_counter,3)
    #train_trigram_witten_bell = train_char_ngram(fourgram_letter_counter,3,fourgram_witten_bell)
    #test_trigram_witten_bell = test_part(fourgram_letter_counter,3)
    #train_trigram_backoff = train_char_ngram(fourgram_letter_counter,3,four_backoff_prob)
    #test_trigram_backoff = test_part(fourgram_letter_counter,3)
    #print sorted(train_bigram_witten_bell.iteritems(),key=lambda (k,v):v,reverse=True)
    #print sorted(test_bigram_add_one.iteritems(),key=lambda (k,v):v,reverse=True)
    
    #ranked_train = convert_train_to_rank(train_unigram_unsmoothing)
    #ranked_train = convert_train_to_rank(train_unigram_add_one)
    #ranked_train = convert_train_to_rank(train_unigram_witten_bell)
    #ranked_train = convert_train_to_rank(train_unigram_backoff)
    #ranked_train = convert_train_to_rank(train_bigram_unsmoothing)
    #ranked_train = convert_train_to_rank(train_bigram_add_one)
    #ranked_train = convert_train_to_rank(train_bigram_witten_bell)
    #ranked_train = convert_train_to_rank(train_bigram_backoff)
    #ranked_train = convert_train_to_rank(train_trigram_unsmoothing)
    #ranked_train = convert_train_to_rank(train_trigram_add_one)
    #ranked_train = convert_train_to_rank(train_trigram_witten_bell)
    #ranked_train = convert_train_to_rank(train_trigram_backoff)
    #print sorted(ranked_train.iteritems(),key=lambda (k,v):v,reverse=True)

    #all_train_bigrams_nr = sum(bigram_letter_counter.values())
    #all_train_trigram_nr = sum(trigram_letter_counter.values())
    #all_train_fourgrams_nr = sum(fourgram_letter_counter.values())
    
    #bigram_test_letter_counter = sum(bigram_letter_counter.values())
    #trigram_test_letter_counter = sum(trigram_letter_counter.values())
    #fourgram_test_letter_counter = sum(fourgram_letter_counter.values())
    
    #result= result_list(test_unigram_unsmoothing,ranked_train,all_train_bigrams_nr)
    #result= result_list(test_unigram_add_one,ranked_train,all_train_bigrams_nr)
    #result= result_list(test_unigram_witten_bell,ranked_train,all_train_bigrams_nr)
    #result= result_list(test_unigram_backoff,ranked_train,all_train_bigrams_nr)
    #result= result_list(test_bigram_unsmoothing,ranked_train,all_train_trigram_nr)
    #result= result_list(test_bigram_add_one,ranked_train,all_train_trigram_nr)
    #result= result_list(test_bigram_witten_bell,ranked_train,all_train_trigram_nr)
    #result= result_list(test_bigram_backoff,ranked_train,all_train_trigram_nr)
    #result= result_list(test_trigram_unsmoothing,ranked_train,all_train_fourgrams_nr)
    #result= result_list(test_trigram_add_one,ranked_train,all_train_fourgrams_nr)
    #result= result_list(test_trigram_witten_bell,ranked_train,all_train_fourgrams_nr)
    #result= result_list(test_trigram_backoff,ranked_train,all_train_fourgrams_nr)

    #print '~' * 80
    #print 'trigram backoff probability'
    #print 'train'+str(i)+' character number: ' + str(all_train_trigram_nr)
    #print 'test'+str(i)+' character number: ' + str(trigram_test_letter_counter)
    #for key,value in result.items():
    #    print 'top' +str(key)+': ' + str(value)+'%'
    
    #for ngram,value in bigram_letter_counter.iteritems():
    #    print 'bigram: "'+str(ngram) +'"'
    #    print 'unsmoothing prob: '+str(bigram_unsmoothing_prob[ngram])
    #    print 'add-one prob: '+str(bigram_add_one_prob[ngram])
    #    print 'witten-bell prob: '+str(bigram_witten_bell[ngram])
    #    print 'backoff prob: '+str(bigram_backoff_prob[ngram])
    #    print '~'*80