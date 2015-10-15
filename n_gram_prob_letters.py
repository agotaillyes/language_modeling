import profile
import itertools
import collections
import string
import re
import sys

# megkapjuk a beolvasott file osszes token-et (szavat) kozpontozas nelkul
def get_tokens_list(file_name):
    tokens_list = []
    punct = set(string.punctuation)
    '!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890'

    with open(file_name,'r') as f:
<<<<<<< HEAD
         for line in f:
            for word in line.split():
                if not word.startswith("-"):
                    word = word.lower()
                    word = unicode(word,'utf-8')
                    word = re.sub('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890]','',word)
                    tokens_list.append(word)
=======
     for line in f:
        for word in line.split():
            word = word.lower()
            word = re.sub('[!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~1234567890]','',word)
            word = unicode(word,'utf-8')
            if not word.startswith('-') and word != '':
                tokens_list.append(word)
>>>>>>> some changes with printing
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
<<<<<<< HEAD

=======
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
>>>>>>> some changes with printing
################### MAIN ############################
if __name__ == '__main__':
    file_name = sys.argv[1]
    # beteszi szonak olyanokat, amiket nem kellene
    tokens_list = get_tokens_list(file_name)
    #print sorted(tokens_list,reverse=False)
    all_letters_nr = get_all_letters_nr(tokens_list)
    #print all_letters_nr
    k=11

    ############## UNIGRAM #################
    unigram_letter_counter = n_gram_letter_counter(1,tokens_list)
    letters_type = len(unigram_letter_counter)
    
    sorted_unigram_letter_counter = sorted(unigram_letter_counter.iteritems(),key=lambda (k,v):v, reverse=True)

    unigram_unsmoothing_prob = ngram_unsmoothing_prob(1,unigram_letter_counter,unigram_letter_counter,all_letters_nr)
<<<<<<< HEAD
=======
    #print 'unigram letter counter: ' + str(sorted(unigram_letter_counter.iteritems(), key=lambda (k, v): v, reverse=True))
    
    #print 'unigram letter counter: ' + str(sorted(changeSmallNrOfCharacters(unigram_letter_counter).iteritems(),key=lambda (k,v): v,reverse=True))
>>>>>>> some changes with printing
    #print 'unigram unsmoothing probabilities: ' + str(sorted(unigram_unsmoothing_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    unigram_add_one_prob = ngram_add_one_prob(1,unigram_letter_counter,unigram_letter_counter,all_letters_nr,letters_type)
    #print 'unigram add-one probabilities: ' + str(sorted(unigram_add_one_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    unigram_witten_bell = ngram_witten_bell_prob(1,unigram_letter_counter,1)
    #print 'unigram witten-bell discounting: ' + str(sorted(unigram_witten_bell.iteritems(),key=lambda (k,v): v,reverse=True))
<<<<<<< HEAD

    #unigramokra talan meg nem erdemes kiszamolni a good-turing-discounting-ot
    #unigram_occuring_list = ngram_occuring_list(unigram_letter_counter,k)

    #unigram_prob_tilde_list = prob_tilde(unigram_letter_counter,unigram_letter_counter,all_letters_nr)
=======

    print 'UNIGRAM'
    print '=' * 80
    for value, key in sorted_unigram_letter_counter:
        print 'unigram: ' + value + ':' + str(key)
        print 'unsmoothing prob: ' + str(unigram_unsmoothing_prob[value])
        print 'add-one smoothing prob: ' + str(unigram_add_one_prob[value])
        print 'witten-bell discounting: ' + str(unigram_witten_bell[value])
>>>>>>> some changes with printing

    ############## BIGRAM ##################
    bigram_letter_counter = n_gram_letter_counter(2,tokens_list)
    #print bigram_letter_counter
    #print 'bigram letter counter: ' + str(sorted(bigram_letter_counter.iteritems(), key=lambda (k,v): v,reverse=True))

    sorted_bigram_letter_counter = sorted(bigram_letter_counter.iteritems(),key=lambda (k,v):v, reverse=True)

    bigram_unsmoothing_prob = ngram_unsmoothing_prob(2,bigram_letter_counter,unigram_letter_counter,all_letters_nr)
    #print 'bigram unsmoothing probabilities: ' + str(sorted(bigram_unsmoothing_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    bigram_add_one_prob = ngram_add_one_prob(2,bigram_letter_counter,unigram_letter_counter,all_letters_nr,letters_type)
    #print 'bigram add-one probabilities: ' + str(sorted(bigram_add_one_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    bigram_witten_bell = ngram_witten_bell_prob(2,bigram_letter_counter,unigram_letter_counter)
    #print 'bigram witten-bell discounting: ' + str(sorted(bigram_witten_bell.iteritems(),key=lambda (k,v): v,reverse=True))
<<<<<<< HEAD

    #bigram_occuring_list = ngram_occuring_list(bigram_letter_counter,k)
    #print 'biram occuring list: ' + str(bigram_occuring_list)
    #l = len(bigram_occuring_list)
    #bigram_good_turing_discounting = ngram_good_turing_discounting(bigram_occuring_list,l)
    #print 'bigram good turing discounting: ' + str(bigram_good_turing_discounting)

    #bigram_prob_tilde_list = prob_tilde(bigram_letter_counter,unigram_letter_counter,all_letters_nr)
    #print bigram_prob_tilde_list

    #bigram_alpha = alpha(bigram_prob_tilde_list,unigram_prob_tilde_list,bigram_letter_counter)
    #print alpha(bigram_prob_tilde_list,bigram_list)
=======

    print 'BIGRAM'
    print '=' * 80
    for value, key in sorted_bigram_letter_counter:
        print 'bigram: ' + value + ':' + str(key)
        print 'unsmoothing prob: ' + str(bigram_unsmoothing_prob[value])
        print 'add-one smoothing prob: ' + str(bigram_add_one_prob[value])
        print 'witten-bell discounting: ' + str(bigram_witten_bell[value])
>>>>>>> some changes with printing

    ############## TRIGRAM #################
    trigram_letter_counter = n_gram_letter_counter(3,tokens_list)
    #print 'trigram letter counter' +  str(sorted(trigram_letter_counter.iteritems(),key=lambda (k,v): v, reverse=True))

    sorted_trigram_letter_counter = sorted(trigram_letter_counter.iteritems(),key=lambda (k,v):v, reverse=True)

    trigram_unsmoothing_prob = ngram_unsmoothing_prob(3,trigram_letter_counter,bigram_letter_counter,all_letters_nr)
    #print 'trigram unsmoothing probabilities: ' + str(sorted(trigram_unsmoothing_prob.iteritems(),key=lambda (k,v): v,reverse=True))

    trigram_add_one_prob = ngram_add_one_prob(3,trigram_letter_counter,unigram_letter_counter,all_letters_nr,letters_type)
    #print 'trigram add-one probabilities: ' + str(sorted(trigram_add_one_prob.iteritems(),key=lambda (k,v):v, reverse=True))

    trigram_witten_bell = ngram_witten_bell_prob(3,trigram_letter_counter,bigram_letter_counter)
<<<<<<< HEAD
    #print 'trigram witten-bell discounting' + str(trigram_witten_bell)

    #trigram_occuring_list = ngram_occuring_list(trigram_letter_counter,k)
    #l = len(trigram_occuring_list)
    #print 'trigram occuring list: ' + str(trigram_occuring_list)
    #trigram_good_turing_discounting = ngram_good_turing_discounting(trigram_occuring_list,l)
    #print 'trigram good-turing discounting: ' + str(trigram_good_turing_discounting)

    #trigram_prob_tilde_list = prob_tilde(trigram_letter_counter,bigram_letter_counter,all_letters_nr)
    #trigram_alpha = alpha(trigram_prob_tilde_list,bigram_prob_tilde_list,trigram_letter_counter)
    #trigram_backoff_prob = ngram_backoff(trigram_list,trigram_letter_counter,trigram_prob_tilde_list,bigram_prob_tilde_list,unigram_prob_tilde_list,bigram_letter_counter,bigram_alpha,unigram_alpha)
    #print trigram_list
    #print trigram_list
    #print sorted(trigram_prob_tilde_list.iteritems(), reverse=True)
    #print bigram_prob_tilde_list
    #print 'trigram BACKOFF: ' + str(trigram_backoff_prob)
##############################################################################
    orig_stdout = sys.stdout
    f=file(sys.argv[2], 'w')
    sys.stdout=f
    
    print 'UNSMOOTHING unigram'
    print '=' * 80
    profile.run('print unigram_unsmoothing_prob')
    
    print 'UNSMOOTHING bigram'
    print '=' * 80
    profile.run('print bigram_unsmoothing_prob')
    
    print 'UNSMOOTHING trigram'
    print '=' * 80
    profile.run('print trigram_unsmoothing_prob')
    
    print 'ADD-ONE unigram'
    print '=' * 80
    profile.run('print unigram_add_one_prob ')
    
    print 'ADD-ONE bigram'
    print '=' * 80
    profile.run('print bigram_add_one_prob')
    
    print 'ADD-ONE trigram'
    print '=' * 80
    profile.run('print trigram_add_one_prob')
    
    print 'WITTEN-BELL unigram'
    print '=' * 80
    profile.run('print unigram_witten_bell')
    
    print 'WITTEN-BELL bigram'
    print '=' * 80
    profile.run('print bigram_witten_bell')
    
    print 'WITTEN-BELL trigram'
    print '=' * 80
    profile.run('print trigram_witten_bell')
    
    sys.stdout=orig_stdout
    f.close()
=======
    #print 'trigram witten-bell discounting' + str(sorted(trigram_witten_bell.iteritems(), key=lambda (k,v):v, reverse=True))

    print 'TRIGRAM'
    print '=' * 80
    for value, key in sorted_trigram_letter_counter:
        print 'trigram: ' + value + ':' + str(key)
        print 'unsmoothing prob: ' + str(trigram_unsmoothing_prob[value])
        print 'add-one smoothing prob: ' + str(trigram_add_one_prob[value])
        print 'witten-bell discounting: ' + str(trigram_witten_bell[value])
##############################################################################
    #orig_stdout = sys.stdout
    #f=file(sys.argv[2],'w')
    #sys.stdout=f

    #print 'UNSMOOTHING unigram'
    #print '=' * 80
    #profile.run('print unigram_unsmoothing_prob')

    #print 'ADD-ONE unigram'
    #print '=' * 80
    #profile.run('print unigram_add_one_prob')

    #print 'WITTEN-BELL unigram'
    #print '=' * 80
    #profile.run('print unigram_witten_bell')

    #print 'UNSMOOTHING bigram'
    #print '=' * 80
    #profile.run('print bigram_unsmoothing_prob')

    #print 'ADD-ONE bigram'
    #print '=' * 80
    #profile.run('print bigram_add_one_prob')

    #print 'WITTEN-BELL bigram'
    #print '=' * 80
    #profile.run('print bigram_witten_bell')

    #print 'UNSMOOTHING trigram'
    #print '=' * 80
    #profile.run('print trigram_unsmoothing_prob')

    #print 'ADD-ONE trigram'
    #print '=' * 80
    #profile.run('print trigram_add_one_prob')

    #print 'WITTEN-BELL trigram'
    #print '=' * 80
    #profile.run('print trigram_witten_bell')

    #sys.stdout=orig_stdout
    #f.close()
>>>>>>> some changes with printing
