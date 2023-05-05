package Network;

import java.util.Objects;

/**
 *
 * @author Ozan Ozisik
 */
public class Node {
    private final String name;
    
    public Node(String name){
        this.name=name;
    }
    
    public String getName(){
        return name;
    }
    
    @Override
    public String toString(){
        return getName();
    }
    
    @Override
    public int hashCode() {
        return name.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Node other = (Node) obj;
        if (!Objects.equals(this.name, other.name)) {
            return false;
        }
        return true;
    }

    
    
}
