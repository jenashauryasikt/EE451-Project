#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <chrono>
#include <limits>
#include <stdlib.h>

using namespace std;

int total_ticks_read = 0;
int total_ticks_processed = 0;
int total_orders_executed = 0;
int total_ticks = 0;
const int delay = 10000;

string companies[30] = {"AAPL", "AMZN", "ARISTA", "ATLASSIAN", "BANKNY", "BLKROCK", "BOFA", "CISCO", "CITI", "DELTA", "EXPEDIA", "FACEBOOK", "GOOGLE", "JPMC", "LYFT", "MICROSOFT", "MORGAN", "NETFLIX", "NVIDIA", "ORACLE", "PALOALTO", "PAYPAL", "QUALCOMM", "SALESFORCE", "SLB", "SOUTHWEST", "SPOTIFY", "TESLA", "UBER", "ZOOM"};

class PriceUpdate {
	public:
		string symbol;
		float ask,bid,ask_vol,bid_vol;
		long int timestamp;
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
		}
};

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
		void addUpdatePrice(PriceUpdate tick) {
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
            float gain_sum_temp = 0;
            int gain_count_temp = 0;
            float loss_sum_temp = 0;
            int loss_count_temp = 0;
			if(updates>50) {
                for(int i=0;i<50;i++) {
                    if(gain[updates-i]!=0) {
                        gain_sum_temp += gain[updates-i];
                        gain_count_temp++;
                    } else {
                        loss_sum_temp += loss[updates-i];
                        loss_count_temp++;
                    }
                }
			}
			avg_gain = gain_sum_temp/gain_count_temp;
			avg_loss = loss_sum_temp/loss_count_temp;
			++updates;
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

// struct Portfolio {
//     double capital;
//     double stocks;
// };

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

// Function to calculate RSI with momentum alpha
// double calculateRSI(const std::vector<TickData>& data, size_t index, size_t period) {
double calculateRSI(SymbolData *data) {
    // Implementation of RSI calculation (not provided for brevity)
    // This would involve calculating average gains and losses over a specified period
    // and using those values to compute the RSI.
    // You may choose to use an existing library for technical indicators.
    // if 
    if(data->avg_loss==0) {
        return 100.0;
    } 
    float RS = data->avg_gain/data->avg_loss;
    
    // Placeholder return value; replace with actual calculation
    return 100.0 - (100.0/(1+RS));
}


// Function to calculate momentum alpha
double calculateMomentumAlpha(SymbolData *data, int period) {
    if (data->updates < period) {
        // Not enough data for calculation
        return 0.0;
    }

    // Calculate the rate of change (momentum alpha)
    double priceChange = data->bid[data->updates-1] - data->bid[data->updates-1-period];
    double pastPrice;
    if (priceChange == 0) {
        pastPrice = data->ask[data->updates-1-period];
        priceChange = data->ask[data->updates-1] - pastPrice;       
    }
    else {
        pastPrice = data->bid[data->updates-1-period];
    }
    
    if (pastPrice != 0.0) {
        return priceChange / pastPrice;
    } 
    // Handle the case where past price is zero (to avoid division by zero)
    else {
        return 0.0;
    }
}


// Function to calculate Fibonacci retracement levels
double calculateFibonacciLevel(double high, double low, double level) {
    double range = high - low;
    return high - (range * level);
}


// Function to set target price using Fibonacci retracement levels
float setTargetPriceUsingFibonacci(SymbolData *data, int period) {
    // Assume previous high and low prices for the Fibonacci calculation
    double previousHigh = data->max_bid; // max till most recent 
    double previousLow = data->min_bid; // min till most recent

    // Fibonacci retracement levels (adjust as needed)
    double fibonacciLevel1 = 0.382;
    double fibonacciLevel2 = 0.618;
    int nTerm = period*delay;

    // Calculate Fibonacci levels based on the previous high and low
    double fibLevel1 = calculateFibonacciLevel(previousHigh, previousLow, fibonacciLevel1);
    double fibLevel2 = calculateFibonacciLevel(previousHigh, previousLow, fibonacciLevel2);

    // Set target price based on Fibonacci levels
    float targetPrice;
    if (data->bid[data->updates-1] < fibLevel1) {
        targetPrice = fibLevel1;
    } else if (data->bid[data->updates-1] < fibLevel2) {
        targetPrice = fibLevel2;
    } else {
        // Set a default target price if the current bid is above the Fibonacci levels
        for(int i=0;i<nTerm;i++) {
            double temp = fibLevel1;
            fibLevel1 = fibLevel1 + fibLevel2;
            fibLevel2 = temp;   
        }
        targetPrice = data->bid[data->updates-1] * (1 + data->max_alpha/period);  // Adjust as needed
    }
    return targetPrice;
}

// Function to calculate risk-reward ratio
double calculateRiskRewardRatio(PriceUpdate &price, SymbolData *data, int period) {
    //Target Price
    double target_price = setTargetPriceUsingFibonacci(data, period);

    // Potential Profit
    double potentialProfit = target_price - price.bid;

    // Potential Loss
    double potentialLoss = -price.bid*data->min_alpha/period;

    // Calculate Risk-Reward Ratio
    double riskRewardRatio = (potentialProfit != 0.0) ? std::abs(potentialProfit / potentialLoss) : 0.0;

    return riskRewardRatio;
}


// Function to execute trades based on RSI, momentum alpha, risk-reward ratio, and smart order routing
bool executeTrades(double rsiValue, double riskRewardRatio, double momentumAlpha, SymbolData *data) {
    // Example threshold for triggering a trade based on RSI and momentum alpha
    const double rsiThreshold = 70.0;
    const double momentumAlphaThreshold = 0;
    bool flag = false;

    // Placeholder for current stock symbol; replace with actual implementation
    // std::string currentSymbol = "AAPL";

    // Example trade execution logic
    if (rsiValue > rsiThreshold && momentumAlpha > momentumAlphaThreshold) {
        // Buy signal
        flag = true;
        double baseOrderSize = data->bid_volume[data->updates-1]/4; // Specify a base order size
        double maxOrderSize = data->bid_volume[data->updates-1];  // Specify a maximum order size

        double orderSize = std::min(baseOrderSize * (1 + riskRewardRatio), maxOrderSize);
        // Adjust portfolio and execute order
        usleep(delay);// trade execution delay    int period = 5;

        portfolio->capital -= orderSize*data->bid[data->updates-1]; // Deduct the order cost from capital
        portfolio->stocks[data->symbol] += orderSize; // Assume all-in strategy (buy as many stocks as possible)
        data->quantity[data->updates-1] = orderSize;

        // Print trade information
        // std::cout << "Buy " << orderSize << " worth of " << data.symbol << " stocks.\n";
    } else if (rsiValue < (100-rsiThreshold) && momentumAlpha < -momentumAlphaThreshold) {
        // Sell signal
        flag = true;
        double baseOrderSize = data->ask_volume[data->updates-1]/4; // Specify a base order size
        double maxOrderSize = data->ask_volume[data->updates-1];  // Specify a maximum order size

        double orderSize = std::min(baseOrderSize * (1 + riskRewardRatio), maxOrderSize);
        // Adjust portfolio and execute order
        usleep(delay);// trade execution delay    int period = 5;

        portfolio->capital += orderSize*data->bid[data->updates-1]; // Deduct the order cost from capital
        portfolio->stocks[data->symbol] -= orderSize; // Assume all-in strategy (buy as many stocks as possible)
        data->quantity[data->updates-1] = -orderSize;

        // Print trade information
        // std::cout << "SELL " << orderSize << " worth of " << data.symbol << " stocks.\n";
    } 
    return flag;
    // else {
    //     // Hold (no trade signal)
    //     // Print hold information
    //     std::cout << "Hold position for " << data.symbol << ".\n";
    // }

    // Additional logic for tracking trades, updating portfolio, etc.
    // ...

    // Note: This is a simplified example, and actual trade execution logic may vary based on your specific strategy.
}

int main() {
    // Specify the name of your CSV file
    string filename = "combined-13.11.2023.csv";

    // Open the CSV file
    std::ifstream file(filename);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    // Create empty stock data lists;
    for(int i=0;i<30;i++) {
        allStocksIndex.insert(make_pair(companies[i], i));
        allStocksData.push_back(new SymbolData(companies[i]));
    }

    // Read the file line by line
    string line;
    getline(file, line);
    int start_time = 1117800000;
    auto accu_start_time = chrono::high_resolution_clock::now();
    auto latency_start = chrono::high_resolution_clock::now();
    auto latency_end = chrono::high_resolution_clock::now();
    auto latency_comp = chrono::duration_cast<chrono::microseconds>(latency_end-latency_start); 
    double RSI, alpha, riskRewardRatio;
    int period = 10; // try 1,2,5

    int count = 0;

    double latency_sum = 0;
    int latency_count = 0;
    while (getline(file, line)) {
        // Use a string stream to parse each line
        istringstream iss(line);
        string value;
        vector<std::string> row;
        bool trade;

        // Split the line into individual values
        while (getline(iss, value, ',')) {
            row.push_back(value.c_str());
        }

        // Use the row data
        PriceUpdate price = PriceUpdate(row[4], strtof(row[0].c_str(),NULL), strtof(row[1].c_str(), NULL), strtof(row[2].c_str(), NULL), strtof(row[3].c_str(), NULL), row[5].c_str());
        total_ticks++;

        auto accu_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(accu_time-accu_start_time);
        // cout<<duration.count()<<" "<<price.timestamp-start_time<<endl;
        if(price.timestamp-start_time >= duration.count()) {
            usleep(1000*(price.timestamp-start_time - duration.count()));
            usleep(delay);// stock information delay
            // latency_start = chrono::high_resolution_clock::now();

            if(duration.count()<600000) {
                allStocksData[allStocksIndex[price.symbol]]->addUpdatePrice(price);
	            total_ticks_read++;
            
                alpha = calculateMomentumAlpha(allStocksData[allStocksIndex[price.symbol]], period);
                allStocksData[allStocksIndex[price.symbol]]->updateAlpha(alpha);
            }


            if(duration.count()>=600000) { // 240000
                latency_start = chrono::high_resolution_clock::now();

                allStocksData[allStocksIndex[price.symbol]]->addUpdatePrice(price);
	            total_ticks_read++;
	        
                alpha = calculateMomentumAlpha(allStocksData[allStocksIndex[price.symbol]], period);
                allStocksData[allStocksIndex[price.symbol]]->updateAlpha(alpha);

                RSI = calculateRSI(allStocksData[allStocksIndex[price.symbol]]);

                riskRewardRatio = calculateRiskRewardRatio(price, allStocksData[allStocksIndex[price.symbol]], period);

                
                trade = executeTrades(RSI, riskRewardRatio, alpha, allStocksData[allStocksIndex[price.symbol]]);
                latency_end = chrono::high_resolution_clock::now();
                total_orders_executed++;
                if (trade){
                    latency_comp = chrono::duration_cast<chrono::microseconds>(latency_end-latency_start);
                    latency_sum += latency_comp.count();
                    latency_count++;
                }
            }


        }
        // std::cout<<duration.count()<<" "<<price.timestamp-start_time<<endl;

        count++;
        // if(count%20000==0) 
        // cout<<count<<endl;
        // if(count == 50000) break;
    }

    double latency_avg=0;
    if(latency_count!=0) 
        latency_avg = (double) latency_sum/latency_count;

    for(auto company: companies) {
        // cout<<allStocksData[allStocksIndex[company]]->updates-1<<endl;
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

    // Close the file
    file.close();
    
    std::cout<<"Average latency is "<<latency_avg<<endl;
    std::cout<<"Total ticks "<<total_ticks<<endl;
    std::cout<<"Total ticks read "<<total_ticks_read<<" Total ticks executed "<<total_orders_executed<<endl;
    std::cout<<"Total capital at the end: "<<portfolio->capital<<endl;
    std::cout<<"All stock quantities are zero:"<<endl;
    for(auto company:companies) {
        std::cout<<company<<" : "<<portfolio->stocks[company]<<endl;
    }

    for(auto company:companies) {
        string file_name = "outputs/" + company + "_output.csv";
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
