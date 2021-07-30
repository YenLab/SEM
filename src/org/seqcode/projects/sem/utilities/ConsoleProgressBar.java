package org.seqcode.projects.sem.utilities;

import java.text.DecimalFormat;

public class ConsoleProgressBar {

    private long minimum = 0; // start of progress bar

    private long maximum = 100; // end of progress bar

    private long barLen = 100; // length of progress bar

    private char showChar = '='; // character used to show progress bar

    private DecimalFormat formater = new DecimalFormat("#.##%");

    public ConsoleProgressBar() {
    }

    public ConsoleProgressBar(long minimum, long maximum,
                              long barLen) {
        this(minimum, maximum, barLen, '=');
    }

    public ConsoleProgressBar(long minimum, long maximum,
                              long barLen, char showChar) {
        this.minimum = minimum;
        this.maximum = maximum;
        this.barLen = barLen;
        this.showChar = showChar;
    }

    /**
     * show the progress bar
     * @param value should not be less than start of bar or larger than end of bar
     */
    public void show(long value) {
        if (value < minimum || value > maximum) {
            return;
        }

        reset();
        minimum = value;
        float rate = (float) (minimum*1.0 / maximum);
        long len = (long) (rate * barLen);
        draw(len, rate);
        if (minimum == maximum) {
            afterComplete();
        }
    }

    private void draw(long len, float rate) {
        System.out.print("Progress: ");
        for (int i = 0; i < len; i++) {
            System.out.print(showChar);
        }
        System.out.print(' ');
        System.out.print(format(rate));
    }

    /**
     * shift the cursor back to the start of line
     */
    private void reset() {
        System.out.print('\r'); 
    } 
    
    /**
     * shift to next line if complete
     */
    private void afterComplete() {
        System.out.print('\n');
    }

    private String format(float num) {
        return formater.format(num);
    }

//	    public static void main(String[] args) throws InterruptedException {
//	        ConsoleProgressBar cpb = new ConsoleProgressBar(0, 100, 30, '#');
//	        for (int i = 1; i <= 100; i++) {
//	            cpb.show(i);
//	            Thread.sleep(100);
//	        }
//	    }

}

