import profile
import itertools
import collections
import string
import re
import sys
from collections import *
from scipy.stats import rankdata
from itertools import tee, islice

# megkapjuk a beolvasott file osszes token-et (szavat) kozpontozas nelkul
def get_tokens_list(input_file_name,output_file_name):
    tokens_list = collections.defaultdict(int)
    punct = set(string.punctuation)
    '!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890'

    letters = set('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890]')

    with open(input_file_name,'r') as infile, open(output_file_name,'w') as outfile:
         for line in infile:
            for word in line.split():
                if not word.startswith('-') and word != '':
                    word = re.sub('[!"#$%&\'()+,-./:;<=>?@[\\]^_`{|}~1234567890]','',word)
                    word = word.lower()
                    #word = unicode(word,'utf-8')
                    if word != '':
                        tokens_list[word] += 1
                    outfile.write(word+' ')
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
            alma=n+' '
            korte+=alma
        length=len(korte)
        value2 = korte[:len(korte)-1]
        szilva[value2]=value
    return szilva

def n_gram_word_counter(j,ngram_word_counter,tokens_list):
    ngram_counter=collections.defaultdict(int)
    if j==1:
        ngram_l = [item for item,value in tokens_list.iteritems()]
    else:
        ngram_l = [' '.join(item) for item in itertools.product(tokens_list,repeat=j)]
    
    for ngram in ngram_l:
        ngram_counter[ngram]=0
 
    for ngram,values in ngram_word_counter.iteritems():
        ngram_counter[ngram]=values

    return ngram_counter

################ UNSMOOTHING probabilities ######################

def ngram_unsmoothing_prob(i,ngram_letter_counter,n_1gram_letter_counter,all_word_nr):
    ngram_prob_list = collections.defaultdict(int)
    prob_list=collections.defaultdict(int)

    if i == 1:
        for ngram in ngram_letter_counter:
            ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0) / all_word_nr

    else:
        for ngram in ngram_letter_counter:
            #print ngram
            #print ngram[0:i-1]
            list=ngram.split()
            n_1gram=list[0]
            if n_1gram_letter_counter[n_1gram] != 0:
                #print ngram_letter_counter[ngram]
                #print n_1gram_letter_counter[ngram[0:i-1]]
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0) / n_1gram_letter_counter[n_1gram]
                #print ngram_prob_list[ngram]
            else:
                ngram_prob_list[ngram] = 0.0
            #print '~'*80
    
    for ngram,value in ngram_prob_list.iteritems():
        if value != 0:
            prob_list[ngram]=value
    return prob_list
######################################################################

#################### ADD-ONE probabilities ############################

def ngram_add_one_prob(i,ngram_letter_counter,n_1gram_letter_counter,all_word_nr,word_type):
    ngram_prob_list = collections.defaultdict(int)
    prob_list = collections.defaultdict(int)
    
    if i == 1:
        for ngram in ngram_letter_counter:
            ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+1.0) / (all_word_nr+word_type)
    else:
        for ngram in ngram_letter_counter:
            list = ngram.split()
            n_1gram=list[0]
            #print ngram
            #print ngram[0:i-1]
            #print n_1gram_letter_counter[ngram[0:i-1]]
            #print n_1gram
            #print n_1gram_letter_counter[n_1gram]
            ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+1.0) / (n_1gram_letter_counter[n_1gram]+word_type)
    
    for ngram,value in ngram_prob_list.iteritems():
        if value != 0:
            prob_list[ngram]=value
    return prob_list

########################################################################

############## WITTEN-BELL DISCOUNTING probabilities ####################

def ngram_witten_bell_prob(i,ngram_letter_counter,n_1gram_letter_counter):
    n_1gram_z_list = collections.defaultdict(int)
    n_1gram_watched_types = collections.defaultdict(int)
    tokens_nr = 0
    ngram_prob_list = collections.defaultdict(int)
    prob_list = collections.defaultdict(int)
    z = 0

    if i != 1:
        for ngram in ngram_letter_counter:
            list = ngram.split()
            n_1gram=list[0]
            if ngram_letter_counter[ngram] == 0:
                n_1gram_z_list[n_1gram] += 1

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
        list=ngram.split()
        n_1gram=list[0]
        if i != 1:
            if(ngram_letter_counter[ngram] == 0):
                ngram_prob_list[ngram] = (n_1gram_watched_types[n_1gram]+0.0)/(n_1gram_z_list[n_1gram]*(tokens_nr+n_1gram_watched_types[n_1gram]))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram]+0.0)/(n_1gram_letter_counter[n_1gram]+n_1gram_watched_types[n_1gram])
        else:
            if(ngram_letter_counter[ngram] == 0):
                ngram_prob_list[ngram] = (watched_types_nr + 0.0) / (z * (tokens_nr + watched_types_nr))
            else:
                ngram_prob_list[ngram] = (ngram_letter_counter[ngram] + 0.0) / (tokens_nr + watched_types_nr)

    for ngram,value in ngram_prob_list.iteritems():
        if value != 0:
            prob_list[ngram]=value
    return prob_list

def c_star(ngram_letter_counter):
    z=0
    c_star=collections.defaultdict(int)
    
    for ngram,value in ngram_letter_counter.items():
        if value == 0:
            z += 1
    wached_types_nr=len(ngram_letter_counter) - z
    tokens_nr = sum(ngram_letter_counter.values())
    
    for ngram,value in ngram_letter_counter.items():
        if value == 0:
            c_star[ngram] = wached_types_nr*tokens_nr/((0.0+z)*(tokens_nr+wached_types_nr))
        else:
            c_star[ngram] = value * tokens_nr /(0.0+tokens_nr+wached_types_nr)
    return c_star
            
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
def unigram_prob_tilde(ungram_letter_counter,all_word_nr):
    unigram_prob_tilde_list = collections.defaultdict(int)
    
    for ngram,value in unigram_word_counter.items():
        unigram_prob_tilde_list[ngram]=value/all_word_nr
    
    return unigram_prob_tilde_list

def prob_tilde(ngram_letter_counter,n_1gram_letter_counter,ngram_witten_bell_prob,tokens_nr,all_word_nr,c_star):
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
        list=ngram.split()
        n_1gram=list[0]
        denominator_part = n_1gram_letter_counter[n_1gram]
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
            list=ngram.split()
            n_1gram=list[1]
            sum_ngram_list += prob_plus1tilde_list[ngram]
            sum_n_1gram_list += prob_tilde_list[n_1gram]                
            #print prob_plus1tilde_list[ngram]
            #print prob_tilde_list[ngram[1:]]
    #print sum_ngram_list
    #print sum_n_1gram_list
    
    alpha = (1.0-sum_ngram_list)/(1-sum_n_1gram_list)

    return alpha

def bigram_backoff(bigram_word_counter,unigram_word_counter,prob_tilde,alpha_n_1gram,all_letter_nr):
    backoff_prob_list = collections.defaultdict(int)

    for ngram,value in bigram_word_counter.items():
        if value > 0:
            backoff_prob_list[ngram]=prob_tilde[ngram]
        else:
            list=ngram.split()
            unigram=list[1]
            backoff_prob_list[ngram]=alpha_n_1gram * unigram_word_counter[unigram] / (0.0+all_letter_nr)

    return backoff_prob_list

def trigram_backoff(trigram_letter_counter,bigram_word_counter,unigram_word_counter,prob_tilde,prob_1tilde,bigram_alpha,unigram_alpha,all_leters_nr):
    trigram_backoff=collections.defaultdict(int)
    
    for ngram, value in trigram_letter_counter.items():
        if value > 0:
            list=ngram.split()
            bigram=list[1]+' '+list[2]
            unigram=list[2]
            trigram_backoff[ngram]=prob_tilde[ngram]
        elif bigram_word_counter[bigram] >0 :
            trigram_backoff[ngram]=bigram_alpha*(0.0+prob_1tilde[ngram[1:]])
        else:
            trigram_backoff[ngram]=unigram_alpha*(0.0+ unigram_word_counter[unigram]/all_word_nr)
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
    return 

#################### REPLACE RARLY OCCURED CHARACTERS WITH * ##################
def list_of_small_nr_of_special_char(unigram_word_counter):
    list_char=''
    
    for letter in unigram_word_counter:
        if unigram_word_counter[letter] <= 30:
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

########################### TRAIN #############################
def train_char_ngram(ngram_word_prob,order,orderplus1_ngram_prob):
    lm = defaultdict(Counter)

    for word,value in ngram_word_prob.iteritems():
        #print word
        data = word
        #print '~'*80
        list=word.split()
        history=list[0]
        char=list[1]
        #history,char =data[0:order],data[order]
        #print history
        #print char
        if char != " " and not (history.startswith(" ") or history.endswith(" ")) and value != 0:                
            word1=word
            lm[history][char]=orderplus1_ngram_prob[word]
            #print word1
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

    for word,value in ngram_letter_counter.iteritems():
        data=word
        list = word.split()
        history=list[0]
        char=list[1]
        if char != " " and not (history.startswith(" ") or history.endswith(" ")) and value !=0:
            alma=word
            test_list[history][char] = ngram_letter_counter[alma]
            #print alma
            #print test_list[history][char]
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
    train_file_name_out = sys.argv[2]
    test_file_name_in=sys.argv[3]
    test_file_name_out=sys.argv[4]
    i=sys.argv[5]
        
    train_token_list=get_tokens_list(train_file_name_in,train_file_name_out)
    test_token_list=get_tokens_list(test_file_name_in,test_file_name_out)
    
    with open(train_file_name_out,"r") as myfile:
        data=myfile.read()
    words = re.findall("\w+",data)
    
    tokens_nr = sum(train_token_list.values())
    
    unigram_word_counter = train_token_list 
    all_word_nr = sum(unigram_word_counter.values())
    
    word_type = len(unigram_word_counter)
    
    bigram_counter = Counter(ngram(words,2))
    bigram_count = convert_ngram(bigram_counter)
    bigram_word_counter=n_gram_word_counter(2,bigram_count,train_token_list)
  
    bigram_star = c_star(bigram_word_counter)
    bigram_witten_bell = ngram_witten_bell_prob(2,bigram_word_counter,unigram_word_counter)
    
    unigram_prob_tilde = unigram_prob_tilde(unigram_word_counter,all_word_nr)
    bigram_prob_tilde = prob_tilde(bigram_word_counter,unigram_word_counter,bigram_witten_bell,tokens_nr,all_word_nr,bigram_star)
    
    unigram_alpha = alpha(unigram_prob_tilde,bigram_prob_tilde,unigram_word_counter,bigram_word_counter)

    bigram_backoff_prob=bigram_backoff(bigram_word_counter,unigram_word_counter,bigram_prob_tilde,unigram_alpha,all_word_nr)
 
    train_unigram_backoff = train_char_ngram(bigram_backoff_prob,1,bigram_backoff_prob)    
    test_unigram_backoff = test_part(bigram_word_counter,1) 

    ranked_train = convert_train_to_rank(train_unigram_backoff)

    all_train_bigrams_nr = sum(train_token_list.values())
    bigram_test_letter_counter = sum(test_token_list.values())

    result= result_list(test_unigram_backoff,ranked_train,all_train_bigrams_nr)

    print '~' * 80
    print 'bigram backoff probability'
    print 'train'+str(i)+' character number: ' + str(all_train_bigrams_nr)
    print 'test'+str(i)+' character number: ' + str(bigram_test_letter_counter)
    for key,value in result.items():
        print 'top' +str(key)+': ' + str(value)+'%'
