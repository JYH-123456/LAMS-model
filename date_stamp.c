#include "snow_model.h"
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>


// 计算两个日期时间之间的时间步数
int calculate_time_steps(struct tm start, struct tm end, int hour_interval) {
    time_t start_time = mktime(&start);
    time_t end_time = mktime(&end);
    
    if (start_time == -1 || end_time == -1) {
        return -1;
    }
    
    double diff_seconds = difftime(end_time, start_time);
    double diff_hours = diff_seconds / 3600.0;
    
    // 计算步数，向上取整以确保包含结束时间
    int steps = (int)ceil(diff_hours / hour_interval) + 1;
    
    return steps;
}


// 生成日期时间列表函数
DateList generate_datetime_list(int start_year, int start_month, int start_day, int start_hour,
                                   int end_year, int end_month, int end_day, int end_hour,
                                   int hour_interval) {
    struct tm start_date = {0};
    struct tm end_date = {0};
    
    // 设置起始日期时间
    start_date.tm_year = start_year - 1900;
    start_date.tm_mon = start_month - 1;
    start_date.tm_mday = start_day;
    start_date.tm_hour = start_hour;
    start_date.tm_min = 0;
    start_date.tm_sec = 0;
    // 设置结束日期时间
    end_date.tm_year = end_year - 1900;
    end_date.tm_mon = end_month - 1;
    end_date.tm_mday = end_day;
    end_date.tm_hour = end_hour;
    end_date.tm_min = 59;
    end_date.tm_sec = 59;
    // 计算时间步总数
    int total_steps = calculate_time_steps(start_date, end_date, hour_interval);
    if (total_steps <= 0) {
        DateList empty_list = {NULL, NULL, NULL, NULL, 0};
        return empty_list;
    }
    // 分配内存
    DateList datetime_list;
    datetime_list.years = (int*)malloc(total_steps * sizeof(int));
    datetime_list.months = (int*)malloc(total_steps * sizeof(int));
    datetime_list.days = (int*)malloc(total_steps * sizeof(int));
    datetime_list.hours = (int*)malloc(total_steps * sizeof(int));
    datetime_list.count = total_steps;
    // 初始化当前时间
    struct tm current_date = start_date;
    time_t current_time = mktime(&current_date);
    time_t end_time = mktime(&end_date);
    
    // 填充年月日小时数组
    for (int i = 0; i < total_steps; i++) {
        // 存储当前日期时间
        datetime_list.years[i] = current_date.tm_year + 1900;
        datetime_list.months[i] = current_date.tm_mon + 1;
        datetime_list.days[i] = current_date.tm_mday;
        datetime_list.hours[i] = current_date.tm_hour;
        // 增加指定小时数
        current_date.tm_hour += hour_interval;
        // 标准化时间（处理日期和月份的跨越）
        current_time = mktime(&current_date);
        // 如果超出结束时间，提前终止循环
        if (current_time > end_time) {
            datetime_list.count = i + 1; // 调整实际步数
            break;
        }
    }
    return datetime_list;
}

// 释放日期时间列表内存
void free_datetime_list(DateList* datetime_list) {
    if (datetime_list->years) free(datetime_list->years);
    if (datetime_list->months) free(datetime_list->months);
    if (datetime_list->days) free(datetime_list->days);
    if (datetime_list->hours) free(datetime_list->hours);
    datetime_list->count = 0;
}

// 打印日期时间列表
void print_datetime_list(const DateList* datetime_list) {
    printf("日期时间列表 (%d 个时间点):\n", datetime_list->count);
    for (int i = 0; i < datetime_list->count; i++) {
        printf("%04d-%02d-%02d %02d:00\n", 
               datetime_list->years[i], 
               datetime_list->months[i], 
               datetime_list->days[i],
               datetime_list->hours[i]);
    }
}


// void main(){
//     DateList datetime_list = generate_datetime_list(
//         2022, 3, 1, 0,  // 起始年月日时
//         2022, 3, 2, 23, // 结束年月日时
//         1               // 时间间隔（小时）
//     );
    
//     // 打印日期时间列表
//     print_datetime_list(&datetime_list);
//     free_datetime_list(&datetime_list);
// }




