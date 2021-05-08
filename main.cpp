#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <chrono>
#include <set>

#define u8 uint8_t
#define u16 uint16_t
#define ll long long

using namespace std;
using namespace std::chrono;

typedef complex<double> cd;
const double PI = acos(-1);

int reverse(int num, int lg_n) {
    int res = 0;
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

void fft(vector<cd> & a, bool invert) {
    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n)
        lg_n++;

    for (int i = 0; i < n; i++) {
        if (i < reverse(i, lg_n))
            swap(a[i], a[reverse(i, lg_n)]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd & x : a)
            x /= n;
    }
}

vector<ll> multiply(vector<ll> const& a, vector<ll> const& b) {
    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    int final_size = a.size() + b.size();
    while (n < final_size) {
        n <<= 1;
    }
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }
    fft(fa, true);

    vector<ll> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = round(fa[i].real());
    }
    vector<ll> res = vector<ll>(result.begin(), result.end() - (result.size() - final_size + 1));
    return res;
}

vector<ll> add(vector<ll> const& a, vector<ll> const& b) {
    int res_size = a.size() > b.size() ? a.size() : b.size();
    vector<ll> result(res_size, 0);
    for (int i = 0; i < res_size; ++i) {
        if (i < a.size() && i < b.size()) {
            result.at(i) = a.at(i) + b.at(i);
        } else if (i >= a.size()) {
            result.at(i) = b.at(i);
        } else if (i >= b.size()) {
            result.at(i) = a.at(i);
        }
    }
    return result;
}

vector<ll> multiply_constant(ll constant, vector<ll> const& a) {
    vector<ll> result(a.size(), 0);
    for (int i = 0; i < a.size(); ++i) {
        result.at(i) = constant * a.at(i);
    }
    return result;
}


pair <vector<u8>, vector<u8>> get_even_odd_vectors(vector <u8> u) {
    vector <u8> u_even;
    vector <u8> u_odd;
    for (int i = 0; i < u.size(); ++i) {
        if (i % 2 == 0) {
            u_even.push_back(u.at(i));
        } else {
            u_odd.push_back(u.at(i));
        }
    }
    return {u_even, u_odd};
}

vector <u8> do_xor(vector <u8> first_v, vector <u8> second_v) {
    vector <u8> xor_v(first_v.size());
    for (int i = 0; i < first_v.size(); ++i) {
        xor_v.at(i) = first_v.at(i) ^ second_v.at(i);
    }
    return xor_v;
}

pair <vector<ll>, vector<ll>> calcA(int n, vector <u8> u) {
    if (n == 1) {
        return {{1, 0},
                {0, 1}};
    }
    int i = u.size();
    if (i % 2 == 0) {
        pair <vector<u8>, vector<u8>> even_odd_u = get_even_odd_vectors(u);
        pair <vector<ll>, vector<ll>> f = calcA(n / 2, do_xor(even_odd_u.first, even_odd_u.second));
        pair <vector<ll>, vector<ll>> g = calcA(n / 2, even_odd_u.second);
        return {add(multiply(f.first, g.first), multiply(f.second, g.second)),
                add(multiply(f.first, g.second), multiply(f.second, g.first))};
    } else {
        vector<u8> u_short = vector<u8>(u.begin(), u.begin() + i - 1);
        pair <vector<u8>, vector<u8>> even_odd_u = get_even_odd_vectors(u_short);
        pair <vector<ll>, vector<ll>> f = calcA(n / 2, do_xor(even_odd_u.first, even_odd_u.second));
        pair <vector<ll>, vector<ll>> g = calcA(n / 2, even_odd_u.second);
        if (u.at(u.size() - 1) == 0) {
            return {multiply(f.first, g.first), multiply(f.second, g.second)};
        } else {
            return {multiply(f.first, g.second), multiply(f.second, g.first)};
        }
    }
}

void generateAllBinaryVectors(vector<pair<bool, vector<int> > >& dynamic_constraints, vector<vector<int>>& vec_u,
                              int n, int *arr, int i) {
    if (i == n) {
        vector<int> temp(n, 0);
        for (int j = 0; j < n; ++j) {
            temp.at(j) = arr[j];
        }
        vec_u.push_back(temp);
        return;
    }
    if (!dynamic_constraints.at(i).first) {
        arr[i] = 0;
        generateAllBinaryVectors(dynamic_constraints, vec_u, n, arr, i + 1);

        arr[i] = 1;
        generateAllBinaryVectors(dynamic_constraints, vec_u, n, arr, i + 1);
    } else {

        arr[i] = 0;
        if (dynamic_constraints.at(i).second.empty()) {
            generateAllBinaryVectors(dynamic_constraints, vec_u, n, arr, i + 1);
        } else {
            for (int ii = 0; ii < dynamic_constraints.at(i).second.size(); ++ii) {
                arr[i] = (arr[i] + arr[dynamic_constraints.at(i).second.at(ii)]) % 2;
            }
            generateAllBinaryVectors(dynamic_constraints, vec_u, n, arr, i + 1);
        }
    }
}

vector<ll> computeWEFNaive(int n, int s, vector<vector<int>> vec_u) {
    vector<ll> a(2 * n + 1);

    for (int i = 0; i < vec_u.size(); ++i) {
        vector<u8> u_short = vector<u8>(vec_u.at(i).begin(), vec_u.at(i).begin() + s);
        pair <vector<ll>, vector<ll>> f = calcA(n, u_short);
        if (vec_u.at(i).at(s) == 0) {
            a = add(a, f.first);
        } else {
            a = add(a, f.second);
        }
    }

    return a;
}

int compareTwoMonomials(vector<int> monomial_a, vector<int> monomial_b) {
    int var_num_a = 0;
    int var_num_b = 0;
    for (int i = 0; i < monomial_a.size(); ++i) {
        if (monomial_a.at(i) == 1) {
            var_num_a++;
        }
        if (monomial_b.at(i) == 1) {
            var_num_b++;
        }
    }
    if (abs(var_num_a - var_num_b) > 1) {
        return 0;
    } else if (abs(var_num_a - var_num_b) == 1) {
        int flag = 0;
        int res = 1;
        for (int i = monomial_a.size() - 1; i >= 0; --i) {
            if (monomial_a.at(i) > monomial_b.at(i)) {
                if (flag == 1) {
                    return 0;
                }
                if (flag == 0) {
                    res = 1;
                }
                flag++;
            } else if (monomial_a.at(i) < monomial_b.at(i)) {
                if (flag == 1) {
                    return 0;
                }
                if (flag == 0) {
                    res = -1;
                }
                flag++;
            }
        }
        return res;
    } else {
        int flag = 0;
        int res = 1;
        for (int i = monomial_a.size() - 1; i >= 0; --i) {
            if (monomial_a.at(i) > monomial_b.at(i)) {
                if (flag == 2) {
                    return 0;
                }
                if (flag == 1 && res == 1) {
                    return 0;
                }
                if (flag == 0) {
                    res = 1;
                }
                flag++;
            } else if (monomial_a.at(i) < monomial_b.at(i)) {
                if (flag == 2) {
                    return 0;
                }
                if (flag == 1 && res == -1) {
                    return 0;
                }
                if (flag == 0) {
                    res = -1;
                }
                flag++;
            }
        }
        return res;
    }
}

vector<vector<int>> get_monomials_order(vector<vector<int>> monomials) {
    vector<vector<int>> monomials_order(monomials.size(), vector<int>(monomials.size(), 0));
    for (int i = 0; i < monomials.size(); ++i) {
        for (int j = 0; j < monomials.size(); ++j) {
            if (i == j) {
                monomials_order.at(i).at(j) = 0;
            } else {
                monomials_order.at(i).at(j) = compareTwoMonomials(monomials.at(i), monomials.at(j));
            }
        }
    }
    return monomials_order;
}

vector<ll> computeWEF(int n, int s, vector<vector<int>> vec_u,
                      vector<vector<int>> monomials_order,
                      vector<int> red_indexes) {
    if (red_indexes.size() == 0) {
        vector<u8> u_zero(s, 0);
        pair <vector<ll>, vector<ll>> f = calcA(n, u_zero);
//        cout << "\nf:\n";
//        for (int i = 0; i < f.first.size(); ++i) {
//            if (f.first.at(i) != 0) {
//                cout << i << ": " << f.first.at(i) << ' ';
//            }
//        }
//        cout << "\n\n";
        return f.first;
    } else {
        int f = red_indexes.at(0);

        vector<int> s_set;

        for (int i = 0; i < red_indexes.size(); ++i) {
            if (monomials_order.at(f).at(red_indexes.at(i)) == 1) {
                s_set.push_back(red_indexes.at(i));
            }
        }
        vector<vector<int>> new_vec_u;
        vector<vector<int>> not_used_vec_u;
        for (int i = 0; i < vec_u.size(); ++i) {
            bool all_zeros = true;
            for (int j = 0; j < s_set.size(); ++j) {
                if (vec_u.at(i).at(s_set.at(j)) == 1) {
                    all_zeros = false;
                }
            }
            if (vec_u.at(i).at(f) == 1) {
                if (all_zeros) {
                    new_vec_u.push_back(vec_u.at(i));
                }
            } else {
                not_used_vec_u.push_back(vec_u.at(i));
            }
        }

        vector<int> new_red_indexes;
        for (int i = 0; i < red_indexes.size(); ++i) {
            if (red_indexes.at(i) != f) {
                new_red_indexes.push_back(red_indexes.at(i));
            }
        }

//        if (f == 15) {
//            cout << new_vec_u.size() << '\n';
//            for (int i = 0; i < new_vec_u.size(); ++i) {
//                for (int j = 0; j < new_vec_u.at(i).size(); ++j) {
//
//                    cout << new_vec_u.at(i).at(j) << ' ';
//                }
//
//                cout << '\n';
//            }
//        }
//        cout << "red size: " << red_indexes.size() << '\n';
//        for (int i = 0; i < red_indexes.size(); ++i) {
//            cout << red_indexes.at(i) << ' ';
//        }
//        cout << '\n';
//        cout << '\n';
//        cout << '\n';
//        cout << "set S size: " << s_set.size() << '\n';
//        for (int i = 0; i < s_set.size(); ++i) {
//            cout << s_set.at(i) << ' ';
//        }
//        cout << '\n';

        vector<ll> result_b = computeWEF(n, s, not_used_vec_u,
                                         monomials_order, new_red_indexes);
        vector<ll> result_c = computeWEFNaive(n, s, new_vec_u);



        vector<ll> result = add(result_b, multiply_constant(pow(2, s_set.size()), result_c));
        return result;
    }
}

vector<vector<int>> get_monomials(int n) {
    vector<vector<int>> monomials;
    for (int i = 0; i < (1 << n); ++i) {
        int curRow = i;
        vector<int> b(n, 0);
        int curBinInd = 0;
        while (curRow > 0) {
            b.at(curBinInd) = curRow % 2;
            curBinInd++;
            curRow /= 2;
        }
        vector<int> monomial(n, 0);
        for (int j = 0; j < n; ++j) {
            monomial.at(j) = 1 - b.at(j);
        }
        monomials.push_back(monomial);
    }
    return monomials;
}


int find_the_most_right_one_pos(std::vector<int> &row_vector) {
    int pos = row_vector.size() - 1;
    while (pos > -1) {
        if (row_vector[pos]) {
            return pos;
        }
        --pos;
    }
    return pos;
}


// compute WEF
int main() {
    cin.tie(nullptr);
    ios_base::sync_with_stdio(false);
    freopen("input1.txt", "r", stdin);
//    freopen("output.txt", "w", stdout);

//    int apply_constraints = 0;
//    cin >> apply_constraints;
//    if (apply_constraints) {
    auto start = high_resolution_clock::now();
    int m, n, f_size, last_frozen_pos;

    cin >> m >> f_size >> n >> last_frozen_pos;
    vector<int> frozen_bits(n, 0);
    vector<int> red_indexes;
    for (int i = 0; i < frozen_bits.size(); ++i) {
        cin >> frozen_bits.at(i);
        if (frozen_bits.at(i) == 0) {
            if (i< last_frozen_pos) {
                red_indexes.push_back(i);
            }
        }
    }


    vector<vector<int>> vec_u;
    vector<vector<int>> dynamic_constraints_matrix(f_size, vector<int>(n));

    for (int i = 0; i < f_size; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> dynamic_constraints_matrix.at(i).at(j);
        }
    }

    vector<pair<bool, vector<int> > > dynamic_constraints(n, {false, {}});

    for (int i = 0; i < dynamic_constraints_matrix.size(); ++i) {
        int r_pos = find_the_most_right_one_pos(dynamic_constraints_matrix.at(i));
        dynamic_constraints.at(r_pos).first = true;
        for (int j = 0; j < r_pos; ++j) {
            if (dynamic_constraints_matrix.at(i).at(j) == 1) {
                dynamic_constraints.at(r_pos).second.push_back(j);
            }
        }
    }

//    for (int i = 0; i < dynamic_constraints.size(); ++i) {
//        cout << dynamic_constraints.at(i).first << ' ';
//        for (int j = 0; j < dynamic_constraints.at(i).second.size(); ++j) {
//            cout << dynamic_constraints.at(i).second.at(j) << ' ';
//        }
//        cout << '\n';
//    }
//    cout << '\n';
//    cout << '\n';
    int arr[last_frozen_pos + 1];
    generateAllBinaryVectors(dynamic_constraints, vec_u, last_frozen_pos + 1, arr, 0);

//    for (int i = 0; i < vec_u.size(); ++i) {
//        for (int j = 0; j < vec_u.at(i).size(); ++j) {
//            cout << vec_u.at(i).at(j) << ' ';
//        }
//        cout << '\n';
//    }

//    vector<ll> result_large1 = computeWEFNaive(n, last_frozen_pos, vec_u);
//    vector<ll> result1 = vector<ll>(result_large1.begin(), result_large1.end() - result_large1.size() + n + 1);
//    cout << "\n";
//    for (int i = 0; i < result1.size(); ++i) {
//        cout << result1.at(i) << '\n';
//    }
//    cout << "\n";
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<seconds>(stop - start);
//    cout << '\n';
//    cout << "Time: " << duration.count() << '\n';
    vector<vector<int>> monomials = get_monomials(m);
    vector<vector<int> > monomials_order = get_monomials_order(monomials);
    vector<ll> result_large = computeWEF(n, last_frozen_pos, vec_u, monomials_order,
                                         red_indexes);
//    cout << "monomials:\n";
//    for (int i = 0; i < monomials.size(); ++i) {
////        cout << i << ": ";
//        for (int j = 0; j < monomials.at(i).size(); ++j) {
//            cout << monomials.at(i).at(j) << ' ';
//        }
//        cout << '\n';
//    }
//    cout << '\n';
//    cout << "monomials order:\n";
//    for (int i = 0; i < monomials_order.size(); ++i) {
//        cout << i << ": ";
//        for (int j = 0; j < monomials_order.at(i).size(); ++j) {
//            cout  << monomials_order.at(i).at(j) << ' ';
//        }
//        cout << '\n';
//    }
//    cout << '\n';
    vector<ll> result = vector<ll>(result_large.begin(), result_large.end() - result_large.size() + n + 1);
    cout << "\n";
    for (int i = 0; i < result.size(); ++i) {
        cout << i << ": " << result.at(i) << '\n';
    }
    cout << "\n";
    auto stop1 = high_resolution_clock::now();


    auto duration1 = duration_cast<seconds>(stop1 - start);
    cout << '\n';
    cout << "Time: " << duration1.count() << '\n';
    return 0;
}
