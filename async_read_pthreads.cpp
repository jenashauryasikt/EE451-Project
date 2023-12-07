#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <chrono>
#include <climits>
#include <stdlib.h>
#include <semaphore.h>
#include <pthread.h>
#include <mutex>

#define num_companies 30

using namespace std;

int max_int = INT_MAX;

int total_ticks_read = 0;
int total_orders_executed = 0;
int total_orders_delayed = 0;
int total_ticks = 0;
double latency_sum = 0;
double latency_count = 0;
mutex totalsUpdateMtx;

string companies[num_companies] = {"AAPL", "AMZN", "ARISTA", "ATLASSIAN", "BANKNY", "BLKROCK", "BOFA", "CISCO", "CITI", "DELTA", "EXPEDIA", "FACEBOOK", "GOOGLE", "JPMC", "LYFT", "MICROSOFT", "MORGAN", "NETFLIX", "NVIDIA", "ORACLE", "PALOALTO", "PAYPAL", "QUALCOMM", "SALESFORCE", "SLB", "SOUTHWEST", "SPOTIFY", "TESLA", "UBER", "ZOOM"};

sem_t readPipe;
mutex readPipeMtx;
bool readFileClosed = false;

int start_time = 1117800000; //Reference time
auto accu_start_time = chrono::high_resolution_clock::now();

int n_threads; // n_threads is the threads used for round robin
mutex portfolioMtx;

class PriceUpdate {
	public:
		string symbol;
		float ask,bid,ask_vol,bid_vol;
		long int timestamp;
        int index;
        double rsi,rrr;
        bool rsi_updated, rrr_updated;
		PriceUpdate(string symbol_p, float ask_p, float bid_p, float ask_vol_p, float bid_vol_p, const char * time) {
			symbol = symbol_p;
			ask = ask_p;
			bid = bid_p;
			ask_vol = ask_vol_p;
			bid_vol = bid_vol_p;
			int year,month,day,hour,minute,second,milisecond;
			sscanf(time, "%4d-%2d-%2d %2d:%2d:%2d.%6d", &year, &month, &day, &hour, &minute, &second, &milisecond);
			milisecond /= 1000;
			timestamp = (long int) milisecond + (second*1000) + (minute*60*1000) + (hour*60*60*1000) + (day*24*60*60*1000);		
            index = -1;		
            rsi = 0;
            rrr = 50;
		}
};

class PriceReadQueue {
    public:
        vector<PriceUpdate> queue;
        int produced, consumed;
        PriceReadQueue() {
            produced = -1;
            consumed = -1;
        }
        void enqueue(PriceUpdate priceTick) {
            ++produced;
            queue.push_back(priceTick);
        }
        PriceUpdate dequeue() {
            ++consumed;
            return queue[consumed];
        }
        bool moreToRead() {
            return consumed!=produced;
        }
};

PriceReadQueue *priceReadQueue = new PriceReadQueue();

class SymbolData {
	public:
		string symbol;
		vector<float> ask,bid,ask_volume,bid_volume,gain,loss,quantity;
		int ups, downs, updates;
		float gain_sum,loss_sum, avg_gain, avg_loss, max_alpha, min_alpha, max_bid, min_bid;
		SymbolData() {}
		SymbolData(string symbol_p) {
			symbol = symbol_p;
			ups = 0;
			downs = 0;
			gain_sum = 0;
			loss_sum = 0;
			updates = 0;
			avg_gain = 0;
			avg_loss = 0;
            max_alpha = numeric_limits<float>::min();
            min_alpha = numeric_limits<float>::max();
            max_bid = numeric_limits<float>::min();
            min_bid = numeric_limits<float>::max();
		}
		SymbolData(SymbolData &t) {
			symbol = t.symbol;
			ask = t.ask;
			bid = t.bid;
			ask_volume = t.ask_volume;
			bid_volume = t.bid_volume;
			gain = t.gain;
			loss = t.loss;
            quantity = t.quantity;
			ups = t.ups;
			downs = t.downs;
			updates = t.updates;
			gain_sum = t.gain_sum;
			loss_sum = t.loss_sum;
			avg_gain = t.avg_gain;
			avg_loss = t.avg_loss;
            max_alpha = t.max_alpha;
            min_alpha = t.min_alpha;
            max_bid = t.max_bid;
            min_bid = t.min_bid;
		}
		int addUpdatePrice(PriceUpdate tick) {
			ask.push_back(tick.ask);
			bid.push_back(tick.bid);
            max_bid = max(max_bid, tick.bid);
            min_bid = min(min_bid, tick.bid);
			ask_volume.push_back(tick.ask_vol);
			bid_volume.push_back(tick.bid_vol);
            quantity.push_back(0);
			if(updates == 0) {
				gain.push_back(0);
				loss.push_back(0);	
			} else {
				if(bid[updates]-bid[updates-1]>0) {
					gain.push_back(bid[updates]-bid[updates-1]);
					loss.push_back(0);
					ups++;
					gain_sum += bid[updates]-bid[updates-1];
				} else if (bid[updates]-bid[updates-1]<0) {
					loss.push_back(bid[updates]-bid[updates-1]);
					gain.push_back(0);
					downs++;
					loss_sum += bid[updates]-bid[updates-1];
				} else if (ask[updates]-ask[updates-1]>0) {
					gain.push_back(ask[updates]-ask[updates-1]);
					loss.push_back(0);
					ups++;
					gain_sum += ask[updates]-ask[updates-1];
				} else {
					loss.push_back(ask[updates]-ask[updates-1]);
					gain.push_back(0);
					downs++;
					loss_sum += ask[updates]-ask[updates-1];
				}
			}
			if(updates>50) {
				if(gain[updates-50]!=0) {
					gain_sum -= gain[updates-50];
					ups--;
				} else {
					loss_sum -= loss[updates-50];
					downs--;
				}	
			}
			avg_gain = gain_sum/ups;
			avg_loss = loss_sum/downs;
			return updates++;
		}

        void updateAlpha(float new_alpha) {
            if(new_alpha > max_alpha) 
                max_alpha = new_alpha;
            if(new_alpha < min_alpha)
                min_alpha = new_alpha;
        }
};

map<string, int> allStocksIndex;
vector<SymbolData*> allStocksData;

int ticksToProcess[num_companies];

mutex ticksToProcessMtx[num_companies];

class Portfolio {
    public:
        double capital;
        map<string, double> stocks;
        Portfolio() {
            capital = 0;
            for(auto company: companies) {
                stocks.insert(make_pair(company,0));
            }
        }
};
Portfolio *portfolio = new Portfolio();

double calculateRSI(string symbol) {
    SymbolData *data = allStocksData[allStocksIndex[symbol]];
    if(data->avg_loss==0) {
        return 100.0;
    } 
    float RS = data->avg_gain/data->avg_loss;
    
    return 100.0 - (100.0/(1+RS));
}

double calculateMomentumAlpha(string symbol, int index, int period) {
    SymbolData *data = allStocksData[allStocksIndex[symbol]];
    if (data->updates < period) {
        return 0.0;
    }

    double priceChange = data->bid[index] - data->bid[index-period];
    double pastPrice;
    if (priceChange == 0) {
        pastPrice = data->ask[index-period];
        priceChange = data->ask[index] - pastPrice;       
    }
    else {
        pastPrice = data->bid[index-period];
    }
    
    if (pastPrice != 0.0) {
        return priceChange / pastPrice;
    } 
    else {
        return 0.0;
    }
}


double calculateFibonacciLevel(double high, double low, double level) {
    double range = high - low;
    return high - (range * level);
}


float setTargetPriceUsingFibonacci(string symbol, int index, int period) {
    SymbolData *data = allStocksData[allStocksIndex[symbol]];

    double previousHigh = data->max_bid; // max till most recent 
    double previousLow = data->min_bid; // min till most recent

    double fibonacciLevel1 = 0.382;
    double fibonacciLevel2 = 0.618;

    double fibLevel1 = calculateFibonacciLevel(previousHigh, previousLow, fibonacciLevel1);
    double fibLevel2 = calculateFibonacciLevel(previousHigh, previousLow, fibonacciLevel2);

    float targetPrice;
    if (data->bid[index] < fibLevel1) {
        targetPrice = fibLevel1;
    } else if (data->bid[index] < fibLevel2) {
        targetPrice = fibLevel2;
    } else {
        targetPrice = data->bid[index] * (1 + data->max_alpha/period); 
    }
    return targetPrice;
}

double calculateRiskRewardRatio(string symbol, int index, int period) {
    SymbolData *data = allStocksData[allStocksIndex[symbol]];

    double target_price = setTargetPriceUsingFibonacci(symbol, index, period);

    double potentialProfit = target_price - data->bid[index];

    double potentialLoss = -data->bid[index]*data->min_alpha/period;

    double riskRewardRatio = (potentialProfit != 0.0) ? std::abs(potentialProfit / potentialLoss) : 0.0;

    return riskRewardRatio;
}

struct FuncParam {
    string symbol;
    int index;
    int period;
    double rsiRes;
    double rrrRes;
};

FuncParam rsi_rrr_params[5000]; 
sem_t rsi_func_sem[5000];
sem_t rrr_func_sem[5000];
sem_t rsi_res_sem[5000];
sem_t rrr_res_sem[5000];
bool func_thread_close[5000];

void * rsiFunc(void * args) {
    int thread_id = *((int *) args);
    while(true) {
        sem_wait(&rsi_func_sem[thread_id]);
        rsi_rrr_params[thread_id].rsiRes = calculateRSI(rsi_rrr_params[thread_id].symbol);
        sem_post(&rsi_res_sem[thread_id]);
    }
    pthread_exit(NULL);
}

void * rrrFunc(void * args) {
    int thread_id = *((int *) args);
    while(true) {
        sem_wait(&rrr_func_sem[thread_id]);
        rsi_rrr_params[thread_id].rrrRes = calculateRiskRewardRatio(rsi_rrr_params[thread_id].symbol, rsi_rrr_params[thread_id].index, rsi_rrr_params[thread_id].period);
        sem_post(&rrr_res_sem[thread_id]);
    }
    pthread_exit(NULL);
}

void* processTick(void * args) {
    int thread_id = *((int *) args);
    const int period = 10;
    const double rsiThreshold = 70.0;
    const double momentumAlphaThreshold = 0;
    int index;
    string symbol;
    double momentumAlpha, rsiValue, riskRewardRatio;
    double baseOrderSize, maxOrderSize, orderSize;
    bool check_to_proceed = false;
    int temp_total_ticks_read = 0, temp_total_orders_executed = 0, temp_total_ticks_observed;
    auto latency_start = chrono::high_resolution_clock::now();
    auto latency_end = chrono::high_resolution_clock::now();
    auto latency_comp = chrono::duration_cast<chrono::microseconds>(latency_end-latency_start);
    double temp_latency_sum = 0;
    double temp_latency_count = 0;
    auto accu_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(accu_time-accu_start_time);

    // Threads for RSI and RRR

    pthread_t rsiThread, rrrThread;

    func_thread_close[thread_id] = false;

    sem_init(&rsi_func_sem[thread_id], 0, 0);
    sem_init(&rsi_res_sem[thread_id], 0, 0);
    
    sem_init(&rrr_func_sem[thread_id], 0, 0);
    sem_init(&rrr_res_sem[thread_id], 0, 0);
    
    int rc = pthread_create(&rsiThread, NULL, rsiFunc, (void *)&thread_id);
    if (rc) { cout << "Error:unable to create thread," << rc << endl; exit(-1); }
    
    rc = pthread_create(&rrrThread, NULL, rrrFunc, (void *)&thread_id);
    if (rc) { cout << "Error:unable to create thread," << rc << endl; exit(-1); }

    while(true) {
        sem_wait(&readPipe);
        readPipeMtx.lock();
        if(readFileClosed && !priceReadQueue->moreToRead()) {
            readPipeMtx.unlock();
            break;
        }
        PriceUpdate price =  priceReadQueue->dequeue();
        readPipeMtx.unlock();
        ++temp_total_ticks_read;

        usleep(10000);// stock information delay

        symbol = price.symbol;
        index = price.index;

        accu_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(accu_time-accu_start_time);
        if(duration.count()<600000) {
            momentumAlpha = calculateMomentumAlpha(symbol, index, period);
            allStocksData[allStocksIndex[symbol]]->updateAlpha(momentumAlpha);
            ++temp_total_ticks_observed;
            continue;
        }

        latency_start = chrono::high_resolution_clock::now();

        rsi_rrr_params[thread_id].symbol = symbol;
        sem_post(&rsi_func_sem[thread_id]);
        rsi_rrr_params[thread_id].index = index;
        rsi_rrr_params[thread_id].period = period;
        sem_post(&rrr_func_sem[thread_id]);
        
        momentumAlpha = calculateMomentumAlpha(symbol, index, period);
        allStocksData[allStocksIndex[symbol]]->updateAlpha(momentumAlpha);

        sem_wait(&rsi_res_sem[thread_id]);
        rsiValue = rsi_rrr_params[thread_id].rsiRes;

        if (rsiValue > rsiThreshold && momentumAlpha > momentumAlphaThreshold) {
            // Buy signal
            baseOrderSize = allStocksData[allStocksIndex[symbol]]->bid_volume[index]/4; // Specify a base order size
            maxOrderSize = allStocksData[allStocksIndex[symbol]]->bid_volume[index];  // Specify a maximum order size

            sem_wait(&rrr_res_sem[thread_id]);
            riskRewardRatio = rsi_rrr_params[thread_id].rrrRes;
            
            orderSize = min(baseOrderSize * (1 + riskRewardRatio), maxOrderSize);

            ticksToProcessMtx[allStocksIndex[symbol]].lock();
            check_to_proceed = index==ticksToProcess[allStocksIndex[symbol]];
            ticksToProcessMtx[allStocksIndex[symbol]].unlock();

            if(!check_to_proceed) {
                continue;
            }

            usleep(10000); // Trade execution delay

            portfolioMtx.lock();
            portfolio->capital -= orderSize*allStocksData[allStocksIndex[symbol]]->bid[index]; // Deduct the order cost from capital
            portfolio->stocks[symbol] += orderSize; // Assume all-in strategy (buy as many stocks as possible)
            allStocksData[allStocksIndex[symbol]]->quantity[index] = orderSize;
            portfolioMtx.unlock();

            latency_end = chrono::high_resolution_clock::now();
            latency_comp = chrono::duration_cast<chrono::microseconds>(latency_end-latency_start);
        
            temp_latency_count++;
            temp_latency_sum += latency_comp.count();
        } else if (rsiValue < (100 - rsiThreshold) && momentumAlpha < -momentumAlphaThreshold) {
            // Sell signal
            baseOrderSize = allStocksData[allStocksIndex[symbol]]->ask_volume[index]/4; // Specify a base order size
            maxOrderSize = allStocksData[allStocksIndex[symbol]]->ask_volume[index];  // Specify a maximum order size

            sem_wait(&rrr_res_sem[thread_id]);
            riskRewardRatio = rsi_rrr_params[thread_id].rrrRes;

            orderSize = std::min(baseOrderSize * (1 + riskRewardRatio), maxOrderSize);

            ticksToProcessMtx[allStocksIndex[symbol]].lock();
            check_to_proceed = index==ticksToProcess[allStocksIndex[symbol]];
            ticksToProcessMtx[allStocksIndex[symbol]].unlock();

            if(!check_to_proceed) {
                continue;
            }

            usleep(10000); // Trade execution delay

            portfolioMtx.lock();
            portfolio->capital += orderSize*allStocksData[allStocksIndex[symbol]]->bid[index]; // Deduct the order cost from capital
            portfolio->stocks[symbol] -= orderSize; // Assume all-in strategy (buy as many stocks as possible)
            allStocksData[allStocksIndex[symbol]]->quantity[index] = -orderSize;
            portfolioMtx.unlock();
            
            latency_end = chrono::high_resolution_clock::now();
            latency_comp = chrono::duration_cast<chrono::microseconds>(latency_end-latency_start);
        
            temp_latency_count++;
            temp_latency_sum += latency_comp.count();
        } else {
            sem_wait(&rrr_res_sem[thread_id]);
        }
        ++temp_total_orders_executed;
    }

    func_thread_close[thread_id] = true;

    totalsUpdateMtx.lock();
    total_ticks_read += temp_total_ticks_read;
    total_orders_executed += temp_total_orders_executed;
    total_orders_delayed += (temp_total_ticks_read - temp_total_ticks_observed - temp_total_orders_executed);
    latency_sum += temp_latency_sum;
    latency_count += temp_latency_count;
    totalsUpdateMtx.unlock();
    pthread_exit(NULL);
}

void readFile() {
    string filename = "combined-13.11.2023.csv";

    std::ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        pthread_exit(NULL);
    }

    // Read the file line by line
    string line;
    getline(file, line);

    start_time = 1117800000; 
    accu_start_time = chrono::high_resolution_clock::now();

    auto accu_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(accu_time-accu_start_time);

    while (getline(file, line)) {
        istringstream iss(line);
        string value;
        vector<std::string> row;

        while (getline(iss, value, ',')) {
            row.push_back(value.c_str());
        }

        PriceUpdate price = PriceUpdate(row[4], strtof(row[0].c_str(),NULL), strtof(row[1].c_str(), NULL), strtof(row[2].c_str(), NULL), strtof(row[3].c_str(), NULL), row[5].c_str());
        total_ticks++;

        accu_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(accu_time-accu_start_time);

        if(price.timestamp-start_time >= duration.count()) {
            usleep(1000*(price.timestamp-start_time - duration.count()));
        }

        int index = allStocksData[allStocksIndex[price.symbol]]->addUpdatePrice(price);
        price.index = index;

        ticksToProcessMtx[allStocksIndex[price.symbol]].lock();
        ticksToProcess[allStocksIndex[price.symbol]] = index;
        ticksToProcessMtx[allStocksIndex[price.symbol]].unlock();

        //Add the tick to queue using lock and semaphore
        readPipeMtx.lock();
        priceReadQueue->enqueue(price);
        readPipeMtx.unlock();
        sem_post(&readPipe);
    }

    // Close the file
    file.close();
    readPipeMtx.lock();
    readFileClosed = true;
    readPipeMtx.unlock();
    for(int i=0;i<n_threads;i++)
        sem_post(&readPipe);
}

int main(int argc, char* argv[]) {
    if(argc>1) {
        n_threads = atoi(argv[1]);
    } else {
        n_threads = 1;
    }

    int rc;

    sem_init(&readPipe, 0, 0);

    // Create empty stock data lists;
    for(int i=0;i<num_companies;i++) {
        allStocksIndex.insert(make_pair(companies[i], i));
        allStocksData.push_back(new SymbolData(companies[i]));
        ticksToProcess[i] = -1;
    }    

    //Start the n threads
    pthread_t threads[n_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(int i=0;i<n_threads;i++) {
        rc = pthread_create(&threads[i], &attr, processTick, (void *) &i);
        if (rc) { cout << "Error:unable to create thread," << rc << endl; exit(-1); }
    }

    readFile();
    
    //Join the n threads
    for(int i=0;i<n_threads;i++) {
        rc = pthread_join(threads[i], NULL);
        if (rc) { cout << "Error:unable to join thread," << rc << endl; exit(-1); }
    }

    pthread_attr_destroy(&attr);

    sem_destroy(&readPipe);

    for(auto company: companies) {
        if(portfolio->stocks[company] < 0) {
            //Buy that stock
            portfolio->capital += allStocksData[allStocksIndex[company]]->bid[allStocksData[allStocksIndex[company]]->updates-1]*portfolio->stocks[company];
            portfolio->stocks[company] = 0;
        } else {
            //Sell that stock
            portfolio->capital += allStocksData[allStocksIndex[company]]->ask[allStocksData[allStocksIndex[company]]->updates-1]*portfolio->stocks[company];
            portfolio->stocks[company] = 0;
        }
    }
    std::cout<<"Latency in Microseconds: "<<latency_sum/latency_count<<endl;
    std::cout<<"Total ticks in raw data "<<total_ticks<<endl;
    std::cout<<"Total ticks read "<<total_ticks_read<<" Total orders executed "<<total_orders_executed<<" Total orders delayed "<<total_orders_delayed<<endl;
    std::cout<<"Total capital at the end: "<<portfolio->capital<<endl;
    std::cout<<"All stock quantities are zero:"<<endl;
    for(auto company:companies) {
        std::cout<<company<<" : "<<portfolio->stocks[company]<<endl;
    }

    for(auto company:companies) {
        string file_name = "outputs_" + to_string(n_threads) + "/" + company + "_output.csv";
        std::ofstream filecsv(file_name);
        filecsv << "ASK,BID,QUANTITY\n";
        for(int i=0;i<allStocksData[allStocksIndex[company]]->updates;i++) {
            filecsv << allStocksData[allStocksIndex[company]]->ask[i];
            filecsv << "," << allStocksData[allStocksIndex[company]]->bid[i];
            filecsv << "," << allStocksData[allStocksIndex[company]]->quantity[i] << "\n";
        }
        filecsv.close();
    }
    return 0;
}
