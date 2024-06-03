public class Points {
    //2-dimensional inputs
    int x;    //horizontal
    int y;    //vertical

    public Points(int x, int y){
        this.x = x;
        this.y = y;
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    //Euclidean distance between this point and point p
    public double getDistance(Points p){
        return Math.sqrt(Math.pow(this.x - p.x, 2) + Math.pow(this.y - p.y, 2));
    }
}


