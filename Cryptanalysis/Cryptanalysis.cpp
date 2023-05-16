#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>

using namespace std;

const double englishLetterFrequency[26] = {
0.0817, 0.0149, 0.0278, 0.0425, 0.1270, 0.0223, 0.0202,
0.0609, 0.0697, 0.0015, 0.0077, 0.0403, 0.0241, 0.0675,
0.0751, 0.0193, 0.0009, 0.0599, 0.0633, 0.0906, 0.0276,
0.0098, 0.0236, 0.0015, 0.0197, 0.0007
};



string decode_block_2(string block, vector<int> perm) {
    string result(block.size(), ' ');
    for (int i = 0; i < block.size(); i++) {
        result[i] = block[perm[i]];
    }
    return result;
}
string decode_block(string block, vector<int> perm) {
    string result(block.size(), ' ');
    for (int i = 0; i < block.size(); i++) {
        result[perm[i]] = block[i];
    }
    return result;
}

string decode_text(string text, vector<int> perm, int block_len) {
    string result(text.size(), ' ');
    for (int i = 0; i < text.size(); i += block_len) {
        string block = text.substr(i, block_len);
        string decoded_block = decode_block_2(block, perm);
        for (int j = 0; j < block_len; j++) {
            result[i + j] = decoded_block[j];
        }
    }
    return result;
}

map<string, double> count_bigrams(string str) {
    map<string, double> freq;
    for (int i = 0; i < str.size() - 1; i++) {
        string bigram(2, ' ');
        bigram[0] = str[i];
        bigram[1] = str[i + 1];
        freq[bigram]++;
    }
    double n = str.size() - 1;
    for (auto& entry : freq) {
        entry.second /= n;
    }
    return freq;
}

double distance(map<string, double>& freq1, map<string, double>& freq2) {
    double dist = 0.0;
    for (auto& entry : freq1) {
        string bigram = entry.first;
        double prob1 = entry.second;
        double prob2 = freq2[bigram];
        dist += abs(prob1 - prob2);
    }
    return dist;
}

vector<int> find_key(string ciphertext, int block_len) {
    vector<string> blocks;
    for (int i = 0; i < ciphertext.size(); i += block_len) {
        blocks.push_back(ciphertext.substr(i, block_len));
    }

    vector<int> indices(block_len);
    for (int i = 0; i < block_len; i++) {
        indices[i] = i;
    }
    vector<vector<int>> perms;
    do {
        perms.push_back(indices);
    } while (next_permutation(indices.begin(), indices.end()));

    map<string, double> natural_freq = {
    {"TH", 0.0271}, {"HE", 0.0239}, {"IN", 0.0203}, {"ER", 0.0178},
    {"AN", 0.0161}, {"RE", 0.0141}, {"ES", 0.0132}, {"ON", 0.0132},
    {"ST", 0.0125}, {"NT", 0.0117}, {"EN", 0.0114}, {"AT", 0.0106},
    {"ED", 0.0106}, {"ND", 0.0105}, {"TO", 0.0104}, {"OR", 0.0099},
    {"EA", 0.0096}, {"TI", 0.0094}, {"AR", 0.0085}, {"TE", 0.0083},
    {"NG", 0.0082}, {"AL", 0.0082}, {"IT", 0.0079}, {"AS", 0.0077},
    {"IS", 0.0073}, {"HA", 0.0069}, {"ET", 0.0068}, {"SE", 0.0068},
    {"OU", 0.0067}, {"OF", 0.0066}, {"LE", 0.0061}, {"SA", 0.0058},
    {"VE", 0.0056}, {"RO", 0.0056}, {"RA", 0.0051}, {"RI", 0.0049}
    };

    vector<int> best_perm;
    double min_dist = numeric_limits<double>::max();
    for (auto perm : perms) {
        map<string, double> ciphertext_freq;
        for (auto block : blocks) {
            string decoded_block = decode_block(block, perm);
            map<string, double> block_freq = count_bigrams(decoded_block);
            for (auto& entry : block_freq) {
                string bigram = entry.first;
                double prob = entry.second;
                ciphertext_freq[bigram] += prob;
            }
        }
        for (auto& entry : ciphertext_freq) {
            entry.second /= blocks.size();
        }
        double dist = distance(ciphertext_freq, natural_freq);
        if (dist < min_dist) {
            min_dist = dist;
            best_perm = perm;
        }
    }

    vector<int> key(block_len);
    for (int i = 0; i < block_len; i++) {
        key[best_perm[i]] = i;
    }
    return key;
}


void CalculateLetterFrequency(string text, double letterFrequency[])
{
    int totalLetters = 0;

    for (int i = 0; i < 26; i++)
    {
        letterFrequency[i] = 0;
    }

    for (int i = 0; i < text.length(); i++)
    {
        char c = text[i];

        if (isalpha(c))
        {
            totalLetters++;
            if (isupper(c))
            {
                c = tolower(c);
            }
            letterFrequency[c - 'a']++;
        }
    }

    for (int i = 0; i < 26; i++)
    {
        letterFrequency[i] /= totalLetters;
        letterFrequency[i] *= 100;
    }
}

double CalculateFrequencyDifference(double letterFrequency[])
{
    double difference = 0;

    for (int i = 0; i < 26; i++)
    {
        difference += pow(letterFrequency[i] - englishLetterFrequency[i], 2);
    }

    return difference;
}

void DecryptCaezar()
{
    string cipherText;
    cout << "Введите зашифрованный текст: ";
    cin >> cipherText;
    int prob;
    string probText = "";
    for (int shift = 1; shift <= 25; shift++)
    {
        string plainText = "";
        double letterFrequency[26];

        for (int i = 0; i < cipherText.length(); i++)
        {
            char c = cipherText[i];

            if (isupper(c))
            {
                c = (c - shift - 65 + 26) % 26 + 65;
            }

            else if (islower(c))
            {
                c = (c - shift - 97 + 26) % 26 + 97;
            }

            plainText += c;
        }

        CalculateLetterFrequency(plainText, letterFrequency);

        double difference = CalculateFrequencyDifference(letterFrequency);

        if (shift == 1) {
            prob = difference;
            probText = plainText;
        }
        else {
            if (difference < prob) {
                prob = difference;
                probText = plainText;
            }
        }

        cout << "Сдвиг: " << shift << ", Открытый текст: " << plainText << ", Разница: " << difference << endl;
    }

    cout << "\n\nНаиболее вероятный исходный текст: " << probText << ", Разница: " << prob << endl;
}


void DecryptVigenere()
{
    string cipherText;
    cout << "Введите зашифрованный текст: ";
    cin >> cipherText;

    double threshold = 0.05; 
    double* averageIndexCoincidence = new double[cipherText.length()];

    vector<int> bestKeyLengths;
    for (int i = 2; i <= 20; i++)
    {
        double sumIOC = 0;
        for (int j = 0; j < i; j++)
        {
            string subText = "";
            for (int k = j; k < cipherText.length(); k += i)
            {
                char c = toupper(cipherText[k]);
                if (c >= 'A' && c <= 'Z')
                {
                    subText += c;
                }
            }
            int n = subText.length();

            double IOC = 0;
            for (int l = 0; l < 26; l++)
            {
                double frequency = 0;
                for (int m = 0; m < n; m++)
                {
                    char c = subText[m];
                    if ('A' + l == c)
                    {
                        frequency++;
                    }
                }
                IOC += frequency * (frequency - 1) / (n * (n - 1));
            }
            sumIOC += IOC;
        }
        double avgIOC = sumIOC / i;
        cout << "avgIOC " << i << ": " << avgIOC << endl;
        if (avgIOC - 0.0667 < threshold && avgIOC - 0.0667 > 0)
        {
            bestKeyLengths.push_back(i);
        }
    }
    double letterFrequency[26];

    for (int i = 0; i < bestKeyLengths.size(); i++)
    {
        int keyLength = bestKeyLengths[i];
        string key = "";
        for (int j = 0; j < keyLength; j++)
        {
            string subText = "";
            for (int k = j; k < cipherText.length(); k += keyLength)
            {
                char c = toupper(cipherText[k]);
                if (c >= 'A' && c <= 'Z')
                {
                    subText += c;
                }
            }

            double difference = INFINITY;
            char bestChar;
            for (int l = 0; l < 26; l++)
            {
                char c = l + 'A';
                string attemptKey = "";
                for (int m = 0; m < subText.length(); m++)
                {
                    char c2 = subText[m];
                    if (isupper(c2))
                    {
                        c2 = (c2 - c + 26) % 26 + 'A';
                    }

                    attemptKey += c2;
                }

                CalculateLetterFrequency(attemptKey, letterFrequency);
                double freqDiff = CalculateFrequencyDifference(letterFrequency);

                if (freqDiff < difference)
                {
                    difference = freqDiff;
                    bestChar = c;
                }
            }
            key += bestChar;
        }

        cout << "Ключ для длины " << keyLength << ": " << key << endl;

        string plainText = "";
        for (int j = 0; j < cipherText.length(); j++)
        {
            char c = toupper(cipherText[j]);
            char k = toupper(key[j % keyLength]);

            if (c >= 'A' && c <= 'Z')
            {
                c = (c - k + 26) % 26 + 'A';
            }

            plainText += c;
        }

        cout << "Расшифрованный текст для длины " << keyLength << ": " << plainText << endl;
    }
}

void DecryptSubstitution_cipher() {

    string cipherText;
    cout << "Введите зашифрованный текст: ";
    cin >> cipherText;

    int textLength = cipherText.length();
    vector<int> possibleBlockLengths; 
    for (int i = 2; i <= textLength / 2; i++) {
        if (textLength % i == 0) {
            possibleBlockLengths.push_back(i);
        }
    }

    if (possibleBlockLengths.empty()) {
        cout << "Невозможно подобрать длину блока для данного текста." << endl;
        return;
    }
    
    for (int i = 0; i < possibleBlockLengths.size(); i++) {
        int block_len = possibleBlockLengths[i];

        vector<int> key = find_key(cipherText, block_len);

        cout << "Key: ";
        for (auto k : key) {
            cout << k + 1 << " ";
        }
        cout << endl;

        string plaintext = decode_text(cipherText, key, block_len);
        cout << "Открытый текст: " << plaintext << endl;
    }
}

int gcd(int a, int b) {
    if (b == 0)
        return a;
    else
        return gcd(b, a % b);
}

void DecryptАffine_cipher() {

    string cipherText;
    cout << "Введите зашифрованный текст: ";
    cin >> cipherText;
    string decrypted_text;
    double max_score = 0;
    int best_a, best_b;

    for (int a = 1; a < 26; a++) {
        if (gcd(a, 26) != 1)
            continue;
        for (int b = 0; b < 26; b++) {
            string temp_text;
            double score = 0;

            for (char c : cipherText) {
                if (isalpha(c)) {
                    int num = (a * (c - 'A' - b + 26) % 26);
                    temp_text += char(num + 'A');
                }
                else {
                    temp_text += c;
                }
            }

            int letter_count[26] = { 0 };
            int total_letters = 0;
            for (char c : temp_text) {
                if (isalpha(c)) {
                    letter_count[c - 'A']++;
                    total_letters++;
                }
            }

            for (int i = 0; i < 26; i++) {
                double expected_frequency = englishLetterFrequency[i];
                double actual_frequency = (double)letter_count[i] / total_letters;
                score += abs(expected_frequency - actual_frequency);
            }

            if (max_score == 0 || score < max_score) {
                max_score = score;
                decrypted_text = temp_text;
                best_a = a;
                best_b = b;
            }
        }
    }

    cout << "Расшифрованный текст: " << decrypted_text << endl;
    cout << "Лучшее значение ключа 1: " << best_a << endl;
    cout << "Лучшее значение ключа 2: " << best_b << endl;
}

int main()
{
    setlocale(LC_ALL, "Rus");
    int option;
    cout << "Выберите метод шифрования:\n1. Шифр Цезаря\n2. Шифр Виженера\n3. Общий шифр перестановки\n4. Аффинный шифр" << endl;
    cin >> option;
    switch (option)
    {
    case 1:
        DecryptCaezar();
        break;
    case 2:
        DecryptVigenere();
        break;
    case 3:
        DecryptSubstitution_cipher();
        break;
    case 4:
        DecryptАffine_cipher();
        break;
    default:
        cout << "Неверный вариант." << endl;
        break;
    }
    return 0;
}